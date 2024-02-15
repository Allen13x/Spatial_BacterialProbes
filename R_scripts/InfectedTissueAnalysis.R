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





# Infected Slice analysis -------------------------------------------------

## Get mabs spot

DefaultAssay(s4data)<-'Spatial'


sum_counts<-s4data@assays$Spatial@counts['rpoB', ] + s4data@assays$Spatial@counts['MAB-0545', ]

s4data$MABS<-sum_counts>0
s4data$MABS_count<-sum_counts


# Subset only positives ---------------------------------------------------






s4p_data<-subset(s4data,MABS & orig.ident!='B1')

DefaultAssay(s4p_data)<-'integrated'





preprocessData(s4p_data,by='integrated')->s4p_data

s4p_data$seurat_clusters<-s4p_data$integrated_snn_res.0.3

Idents(s4p_data)<-'seurat_clusters'


SpatialDimPlot(s4p_data,group.by = 'integrated_snn_res.1',crop=F,combine = F,pt.size.factor = 1)->p


wrap_plots(p[[3]],p[[4]])

s4p_data$fcluster<-factor(s4p_data$seurat_clusters,levels = c('0','2','1','3'))
levels(s4p_data$fcluster)<-list('Upper airways'='0',
                                'Upper airways Submucosal'='2',
                                'Lower airways Parenchima'='1',
                                'Lower airways Granuloma-like'='3')

SpatialDimPlot(s4p_data,group.by = 'fcluster',crop=F,combine = F,pt.size.factor = 1.2,stroke = 0.2,image.alpha = 0.65)->p


wrap_plots(p[[3]]+
             scale_fill_manual(values=c('yellow','purple','#60b038','red')),
           
           p[[4]]+
             scale_fill_manual(values=c('yellow','purple','#60b038','red')))


# Distances Plot ----------------------------------------------------------

DefaultAssay(s4data)='RCTD'

s<-c(A1='A1',B1='B1',C1='C1',D1='D1')
s[3:4]
lapply(s[3:4],function(x){
  
  s4data@images[[x]]@coordinates %>% 
    rownames_to_column('spotID') %>% 
    left_join(s4data@meta.data %>%dplyr::select(MABS) %>%  rownames_to_column('spotID'))->a1
  
  a1 %>% 
    filter(MABS) %>% 
    select(row1=imagerow,col1=imagecol,mabsSpot=spotID) %>% 
    mutate(a='a')->b1
  
  
  a1 %>% 
    mutate(a='a') %>% 
    left_join(b1) %>% 
    mutate(dist=sqrt((row1-imagerow)^2+(col1-imagecol)^2))->c1
  return(c1)
}) %>% bind_rows(.id='Slice') %>%
  select(spotID,dist,Slice,mabsSpot) %>% 
  left_join(FetchData(s4data,rownames(s4data)) %>% rownames_to_column('spotID')) %>% 
  gather(Type,perc,-spotID,-dist,-Slice,-mabsSpot)->dist_m


cellID=c('B cells',
         'Dendritic cells',
         'Macrophages',
         'Monocytes',
         'Neutrophils',
         'T cells')

lapply(names(table(s4p_data$fcluster)),function(x){
  dist_m %>% 
    filter(dist<1000) %>% 
    mutate(dist=(round(dist,0))) %>% 
    select(spotIDMabs=mabsSpot,dist,perc,Type) %>% 
    group_by(spotIDMabs,dist,Type) %>% 
    summarise(perc=mean(perc)) %>% 
    
    left_join(s4p_data@meta.data %>% rownames_to_column('spotIDMabs') %>% select(spotIDMabs,cluster=fcluster))  %>% 
    filter(Type%in%cellID) %>% 
    filter(!Type %in%c('NK cells','Erythrocytes')) %>% 
    filter(cluster==x) %>%
    mutate(cluster=str_wrap(cluster,20)) %>% 
    ggplot(aes(dist, perc, color=Type)) + 
    geom_smooth(method = "loess",lwd=1.5) + scale_color_manual(values=c(gg_color_hue(11)[c(2,4)],'red',gg_color_hue(11)[c(8,9,11)])) + 
    facet_wrap(~cluster)+
    theme_classic()+
    theme(strip.background=element_rect(colour="transparent",
                                        fill="transparent"),
          strip.text = element_text(size=20),
          axis.text=element_text(size=15),
          legend.title = element_blank(),
          axis.title=element_blank())+
    theme(legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          #axis.title = element_text(size=15),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent",
                                         color = NA))+
    guides(fill=F)+
    scale_y_continuous(labels = scales::percent)+
    coord_cartesian(ylim=c(0, 0.41))+
    ylab('Proportion of cell types')+xlab(expression(paste('Distance from M.abs probes(',mu,'m)')))
})->p

ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],ncol=2,nrow=2,common.legend = T,legend='top')->f


annotate_figure(f,left=text_grob("Proportion of cell types", rot = 90,size=20),
                bottom=text_grob(size=20,expression(paste('Distance from MABS+ (',mu,'m)'))))&
  theme(legend.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        axis.title = element_text(size=30),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       color = NA))

# MABS+ vs MABS- -------------------------------------------------------------------------


DefaultAssay(s4data)<-'SCT'

s4data$degclass=case_when(
  (colnames(s4data)%in%n_mabs & s4data$Infection!='CTRL')~1,
  (!colnames(s4data)%in%n_mabs & s4data$Infection!='CTRL')~0,
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



# Probe UMI counts by zone ------------------------------------------------

s4p_data$fcluster=s4p_data$seurat_clusters
data.frame(
  orig.ident=s4p_data$orig.ident,
  cluster=s4p_data$fcluster,
  rpob=s4p_data@assays$Spatial@counts[c('rpoB'),],
  MAB0545=s4p_data@assays$Spatial@counts[c('MAB-0545'),]) %>%
  mutate(orig.ident=ifelse(orig.ident=='A1','A1','merged')) %>%
  filter(rpob+MAB0545>0) %>% group_by(orig.ident,cluster) %>%
  gather(gene,count,-cluster,-orig.ident) %>% group_by(orig.ident,cluster,gene) %>% 
  summarise(c=sum(count)) %>% mutate(cc=sum(c)) %>%distinct() %>% 
  ungroup() %>% 
  mutate(gene=factor(gene,levels=c('rpob','MAB0545'))) %>%
  ggplot()+
  geom_col(aes(x=cluster,y=c,fill=gene),position = 'dodge',col='black',width=0.8) +

  scale_fill_manual(values=rev(c('#D95F02','#1B9E77')),name='',labels=rev(c('lsr2','rpoB')))+
  theme_classic()+
  scale_x_discrete(labels=str_wrap(c('Upper Airways','Upper Airways - Submucosal','Lower Airways - Parenchima','Lower Airways - Granuloma-like structure'),15))+
  ylab('UMI counts')+
  theme(axis.text.x = element_text(angle=45,vjust=0.6),
        axis.text = element_text(size=15),
        axis.title.y = element_text(size=20),
        axis.title.x = element_blank())

# Save Objects -------------------------------------------------------------



saveRDS(s4data,'Path/To/saved/object.rds')
saveRDS(s4p_data,'Path/To/saved/subsetted/object.rds')
