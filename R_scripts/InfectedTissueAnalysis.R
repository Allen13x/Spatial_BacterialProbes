# Load Libraries ----------------------------------------------------------

library(Seurat)
library(patchwork)
library(rstatix)
library(ggpubr)
library(clusterProfiler)
library(RColorBrewer)
library(fgsea)
library(pheatmap)
library(tidyverse)



# Define functions --------------------------------------------------------



preprocessData<-function(x_data,by='SCT'){
  
  DefaultAssay(x_data)<-by
  
  x_data <- RunPCA(x_data, assay = by, verbose = FALSE)
  x_data <- FindNeighbors(x_data, reduction = "pca", dims = 1:30)
  x_data <- RunUMAP(x_data, reduction = "pca", dims = 1:30)
  
  for (i in seq(0.1,1,0.1)){
    x_data<-FindClusters(x_data,resolution = i)
  }
  
  print(SpatialDimPlot(x_data,group.by = paste(by,'_snn_res.',seq(0.1,1,0.1),sep=''),combine = F))
  return(x_data)
  
}



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}




# Load Object -------------------------------------------------------------

s4data<-readRDS('Path/To/saved/object.rds')


# Clusters palette --------------------------------------------------------

Cluster_palette=c(
  
  '#d67d0f',
  '#6bd408',
  '#0897d4',
  '#decc21',
  '#e60b0b',
  '#731024',
  '#B08CA3',
  'grey')


s=levels(s4data$Clusters)
names(Cluster_palette)=s



# Infected Slice analysis -------------------------------------------------

## Get mabs spot

DefaultAssay(s4data)<-'Spatial'


sum_counts<-s4data@assays$Spatial@counts['rpoB', ] + s4data@assays$Spatial@counts['MAB-0545', ]

s4data$MABS<-sum_counts>0
s4data$MABS_count<-sum_counts


# Inflammatory clusters cell composition ----------------------------------

DefaultAssay(s4data)<-'RCTD'  
s4data@meta.data %>% 
  rownames_to_column('spot') %>%  
  left_join(
    FetchData(s4data,rownames(s4data)) %>% 
      rownames_to_column('spot') %>%
      gather(cell,perc,-spot)) %>% 
  
  select(spot,var=matches('cell'),
         splitvar=matches('Clusters',ignore.case=F),Infection,perc) %>%
  mutate(var=as.factor(var))->df


df %>% 
  group_by(splitvar,var,Infection) %>% 
  filter(splitvar%in%c('I1','I2','I3')) %>% 
  mutate(perc=filter_lims(perc)) %>% 
  filter(!var%in%c('NK cells' ,'Erythrocytes','Club cells','Endothelial cells','Fibroblast')) %>% 
  mutate(var=str_wrap(var,15)) %>% 
  ggplot(aes(x=var,y=perc)) +
  geom_boxplot(aes(fill=splitvar),na.rm = T,coef=7,lwd=0.5,show.legend = T)+
  coord_flip()+
  theme_classic()+
  scale_fill_manual(values=Cluster_palette)+
  ggtitle(str_to_title(str_replace_all(var,'_',' ')))+
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=10),
        axis.text.x =element_text(size=10), 
        axis.title = element_blank())+
  scale_y_continuous(labels = scales::percent)





# MABS+ vs MABS- -------------------------------------------------------------------------


DefaultAssay(s4data)<-'SCT'

s4data$degclass=case_when(
  (s4data$MABS & s4data$Infection!='CTRL')~1,
  (s4data$MABS & s4data$Infection!='CTRL')~0,
  s4data$Infection=='CTRL'~2)

Idents(s4data)<-'degclass'
FindMarkers(s4data,ident.1 = 1,ident.2 = 0)->NMABS_markers





# Bubble Plot -------------------------------------------------------------



DotPlot(s4data,idents=c(0,1),features =  c(NMABS_markers %>%
                                             rownames_to_column('gene') %>% 
                                             slice_max(avg_log2FC,n=10) %>%
                                             filter(gene!='Gp2') %>% 
                                             pull(gene) %>% unique(),'S100a9',
                                           NMABS_markers %>%
                                             rownames_to_column('gene') %>% 
                                             slice_min(avg_log2FC,n=10) %>% pull(gene) %>% unique()))+
  theme(axis.text.x = element_text(angle=45,hjust=1),
        axis.title=element_blank())+
  scale_color_gradient(low='grey',high='red')+
  scale_y_discrete(label=c('MABS -','MABS +'))+
  theme(legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       color = NA))


# Enriched cellular populations by cluster --------------------------------

DefaultAssay(s4data)<-'RCTD'

s4data@meta.data %>% 
  rownames_to_column('spotID') %>% 
  filter(Clusters!='C1') %>% 
  filter(Infection!='Ctrl') %>% 
  tibble() %>% 
  mutate(lsr2=lsr2>0 & Infection!='Ctrl',
         lsr2p=mabs & lsr2>0 & Infection!='Ctrl',
         lsr2n=mabs & lsr2<1 & Infection!='Ctrl',
         rpoB=rpoB>0 & Infection!='Ctrl',
         mabsn=!mabs) %>% 
  dplyr::select(spotID,Clusters,orig.ident,mabs) %>%
  left_join(FetchData(s4data,rownames(s4data)) %>% 
              rownames_to_column('spotID')) %>% 
  gather(cell,p,-spotID,-mabs,-Clusters,-orig.ident) %>%
  mutate(mabs=ifelse(mabs,'+','-')) %>% 
  group_by(orig.ident,Clusters,cell) %>% 
  wilcox_test(p~mabs,detailed=T)%>% 
  ungroup() %>% adjust_pvalue(method = 'holm') ->dbs


dbs %>%   
  group_by(Clusters,cell) %>% 
  summarise(agree=prod(estimate)>0,Estimate=str_c(paste(orig.ident,estimate,sep='_'),collapse=';'),p=str_c(p,collapse=',')) %>% 
  filter(agree) %>% 
  group_by(Clusters,cell) %>% 
  separate(p,c('p1','p2'),sep=',') %>% 
  mutate(fisher=MADAM::fisher.method(data.frame(c(as.numeric(p1)[1]),c(as.numeric(p2)[1])))$p.value) %>% 
  ungroup() %>% 
  adjust_pvalue(p.col = 'fisher',method='fdr') %>% 
  filter(fisher<0.05) %>% 
  extract(Estimate,c('d','C1','dd','D1'),'(.*)_(.*);(.*)_(.*)') %>% 
  mutate(Mabs_up=ifelse(C1<0,-log(fisher.adj),log(fisher.adj)),
         Mabs_up=ifelse(Mabs_up==Inf,16,Mabs_up),
         MD=(as.numeric(C1)+as.numeric(D1))/2) %>% 
  dplyr::select(-d,-D1,-dd,-C1) %>% 
  ungroup() %>% 
  dplyr::select(Clusters,cell,Mabs_up,MD)->DD


s4data@meta.data %>% 
  rownames_to_column('spotID') %>% 
  filter(Infection!='Ctrl') %>% 
  dplyr::select(spotID,Clusters,orig.ident,mabs) %>% 
  left_join(
    FetchData(s4data[,s4data$Infection!='Ctrl'],vars = rownames(s4data)) %>% 
      rownames_to_column('spotID')) %>%
  gather(cell,p,-spotID,-Clusters,-orig.ident,-mabs) %>% 
  group_by(Clusters,cell) %>% 
  summarise(p=median(p))->DC

DD %>%  
  left_join(DC) %>% 
  filter(Clusters!='I2') %>% 
  filter(cell%in%c('B cells','Macrophages','Monocytes','Neutrophils','Dendritic cells','T cells')) %>% 
  ggplot()+
  geom_point(aes(x=factor(cell,levels=c('Macrophages','Monocytes','B cells','Neutrophils','Dendritic cells','T cells')),y=factor(Clusters,levels=rev(c('I1','I3','S1','S2','S3'))),size=p,col=Mabs_up<0))+
  theme_classic()+
  xlab('Cells')+
  ylab('Clusters')+
  scale_colour_manual(name='Enriched in',labels=c('MABS +','MABS -'),values=(c('red','royalblue')))+
  scale_size(name='Median\nCellular Proportion')+
  
  guides(colour=guide_legend(nrow=2,override.aes=list(size=3,shape=16)),
         size=guide_legend(nrow=2,byrow=T))+
  theme(axis.text.x = element_text(angle=45,hjust=1,size=10),
        legend.position = 'top',
        plot.margin = margin(10, 100, 20, 100),
        axis.text.y=element_text(size=12))


# Subset only positive spots ----------------------------------------------

s4p_data=subset(s4data,MABS)


# Probe UMI counts --------------------------------------------------------

data.frame(
  orig.ident=s4p_data$orig.ident,
  rpob=s4p_data@assays$Spatial@counts[c('rpoB'),],
  MAB0545=s4p_data@assays$Spatial@counts[c('MAB-0545'),]) %>%
  mutate(orig.ident=ifelse(orig.ident=='A1','A1','merged')) %>%
  filter(rpob+MAB0545>0) %>% 
  group_by(orig.ident) %>% add_count() %>% 
  gather(gene,count,-n,-orig.ident) %>% group_by(orig.ident,gene) %>% 
  summarise(c=sum(count)) %>% mutate(cc=sum(c)) %>%distinct() %>% 
  ungroup() %>% 
  mutate(gene=factor(gene,levels=c('rpob','MAB0545'))) %>% 
  ggplot()+
  geom_col(aes(x=gene,y=c,fill=gene),col='grey11',width=0.5) +

  scale_fill_manual(values=rev(c('#D95F02','#1B9E77')),name='',labels=rev(c('lsr2','rpoB')))+

  theme_classic()+
  ylab('UMI counts')+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.text = element_text(size=15),
    axis.title.y = element_text(size=20),
    axis.title.x = element_blank())+
  theme(legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       color = NA))+
  theme(legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       color = NA))





# INOS Violin plot --------------------------------------------------------

library(gridExtra)
library(grid)




DefaultAssay(s4p_data)<-'SCT'

data.frame(
  orig.ident=s4data$orig.ident,
  rpob=s4p_data@assays$Spatial@counts[c('rpoB'),],
  MAB0545=s4p_data@assays$Spatial@counts[c('MAB-0545'),]) %>% 
  mutate(mabs=rpob+MAB0545>0) %>% 
  rownames_to_column('ID') %>% 
  left_join(data.frame(cluster=s4p_data$fcluster,
                       INOS=s4p_data$IPOSSIA1,
                       orig.ident=s4p_data$orig.ident,
                       spotID=s4p_data$near_mabs) %>% 
              rownames_to_column('ID')) %>% 
  filter(spotID==1) %>% 
  gather(gene,count,rpob,MAB0545) %>%
  mutate(cc=ifelse(count>0,1,0)) %>% 
  filter(count!=0) %>%
  group_by(ID) %>% 
  mutate(gene=ifelse(n()>1,'MAB0545',gene)) %>% 
  ungroup() %>% 
  mutate(gene=factor(gene,levels=c('rpob','MAB0545'))) %>% 
  ungroup() %>% 
  distinct() ->dvin


dvin %>% 
  ggplot()+
  geom_violin(aes(x=gene,y=INOS,fill=gene),shape=21)+
  geom_boxplot(aes(x=gene,y=INOS),width=0.4,shape=21)+
  theme_classic()+
  ylab('INOS pathway levels')+
  scale_fill_manual(values=c('royalblue','#D95F02'),name='',labels=c('MABS+ lsr2-','MABS+ lsr2+'))+
  theme(axis.text.x=element_blank(),
        legend.position = 'None',
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank())->dv


data.frame(lsr2=c('+','-'),MABS=c('+','+')) %>% 
  gather(gene,sign) %>% 
  mutate(y=ifelse(gene=='MABS',2,1)) %>% 
  mutate(y=c('MABS','lsr2','lsr2','MABS')) %>% 
  ggplot()+
  geom_text(aes(y=y,x=gene,label=sign),size=10)+
  geom_vline(xintercept = 1.5,lty=2)+
  theme_void()+
  theme(axis.text.y = element_text(size=15))->t
arrangeGrob(dv,t,heights=c(3,.5))->g



# Save Objects -------------------------------------------------------------



saveRDS(s4data,'Path/To/saved/object.rds')
saveRDS(s4p_data,'Path/To/saved/subsetted/object.rds')
