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
s4p_data<-readRDS('Path/To/saved/subsetted/object.rds')




# Probes map --------------------------------------------------------------



s4data$rpoB<-s4data@assays$Spatial@counts['rpoB', ]
s4data$Lsr2<-s4data@assays$Spatial@counts['MAB-0545', ]


s4data@meta.data %>% 
  rownames_to_column('spot') %>% 
  tibble() %>% 
  mutate(MABS_count=rpoB+Lsr2) %>% 
  column_to_rownames('spot')->s4data@meta.data





DefaultAssay(s4data)<-'Spatial'

p<-list()

for (i in c('rpoB','Lsr2','MABS_count')){
  SpatialFeaturePlot(s4data,features=i,image.alpha = 0.7,crop=F,pt.size.factor = 1.2,stroke=0,combine=F)->p[[i]]}

lapply(c('rpoB','Lsr2','MABS_count'),function(i){
  lapply(p[[i]],function(x){
    x+
      scale_alpha_continuous(range(0,1),limits=c(0,1))+
      scale_fill_gradient(low='dodgerblue',high='blue',limits=c(1,case_when(i=='rpoB'~6,
                                                                            i=='Lsr2'~3,
                                                                            i=='MABS_count'~7)))
  })})->p1


wrap_plots(p1[[1]][[1]],p1[[1]][[2]],p1[[1]][[3]],p1[[1]][[4]],
           p1[[2]][[1]],p1[[2]][[2]],p1[[2]][[3]],p1[[2]][[4]],
           p1[[3]][[1]],p1[[3]][[2]],p1[[3]][[3]],p1[[3]][[4]],ncol=4)



# Define Modules ----------------------------------------------------------


INOS_markers<-str_split('Acp5,Cav1,Cyb5b,Cyb5r3,Cyp1b1,Ddah1,H2-M3,Hsp90aa1,Ifng,Klrk1,Mtarc1,Mtarc2,Nos1,Nos1,Nos1,Nos1,Nos2,Nos2,Nos2,Nos2,Nos2,Nos3,Nos3,Nos3,Nos3,P2rx4,P2rx4,Rora,Slc7a2,Slc7a6,Spr,Ticam1,Tlr2,Tlr2,Tlr4',pattern = ',')[[1]]

ROS_markers<-str_split('Atpif1,Bcl2,Bnip3,Ccn1,Ccn2,Cyb5r4,Cyba,Ddit4,Drd5,Duox1,Gls2,Il19,Il22,Lrrk2,Met,Mpv17,Mpv17l,mt-Nd2,Ndufa13,Ndufs1,Ndufs3,Ndufs4,Nnt,Norad,Nox4,P2rx7,Pax2,Pdgfb,Pdk4,Pla2r1,Pmaip1,Prdx4,Prex1,Rfk,Ripk3,Sesn1,Sesn2,Sh3pxd2a,Snord32a,Snord33,Snord34,Snord35a,Sod1,Tigar,Trp53,Ucp2,Vav1',pattern=',')[[1]]

Hypoxia_markers<-str_split('Acaa2,Adam8,Ado,Adrb2,Ak4,Akt1,Aqp1,Aqp3,Bad,Bcl2,Bmyc,Bnip3,Bnip3l,Cbs,Cpeb1,Cpeb1,Cpeb2,Cpeb2,Cr1l,Egln1,Egln2,Egln3,Eif4ebp1,Epas1,Epas1,Fabp1,Fam162a,Fam162a,Fam162a,Fmn2,Fndc1,Gata6,Gnb1,Gngt1,Hif1a,Hif1a,Hif1a,Hif1a,Hif1a,Hif1a,Hif3a,Higd1a,Higd1a,Higd1a,Hp1bp3,Hyou1,Hyou1,Kcnd2,Kcnk3,Mgarp,Mir199a-2,Mir199a-2,Mir199a-2,Mir214,Mir214,Mir214,Mir668,Mir874,Mlst8,Mtor,Mtor,Myc,Ndnf,Ndp,Ndp,Nkx3-1,Nol3,Nop53,Notch1,Npepps,Oprd1,P4hb,Pgk1,Pink1,Pink1,Ppard,Ppard,Ppard,Pparg,Prkce,Pten,Ptgis,Rgcc,Rora,Rptor,Rtn4,Scn2a,Sdhd,Sirt1,Sirt2,Slc2a4,Slc8a3,Slc9a1,Stub1,Suv39h1,Suv39h2,Tbl2,Tert,Tigar,Trem2,Twist1,Twist1,Ubqln1,Vasn,Vasn,Vegfa,Vegfa,Vldlr,Zfp36l1',pattern = ',')[[1]]


AddModuleScore(s4data,features = list(c(INOS_markers),c(ROS_markers),c(Hypoxia_markers)),name = 'IPOSSIA')->s4data

AddModuleScore(s4p_data,features = list(c(INOS_markers),c(ROS_markers),c(Hypoxia_markers)),name = 'IPOSSIA')->s4p_data


# INOS markers map --------------------------------------------------------

DefaultAssay(s4data)<-'SCT'

SpatialFeaturePlot(s4data,'IPOSSIA1',crop=F,image.alpha = 1,pt.size.factor = 1.2,stroke=0.2,
                   alpha=c(0,1),
                   combine=F)->p
p
lapply(p,function(x){
  x+
    scale_fill_gradient(
      oob=scales::squish,
      limits=c(-0.3,0.2),
      low='lightgreen',high='red',name=str_wrap(str_to_title('INOS pathway levels'),20))+
    theme(legend.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent",
                                         color = NA))->q
  return(q)
})->q

wrap_plots(q[[1]],q[[2]],q[[3]],q[[4]])



# INOS Boxplot Lsr2+ vs Lsr2-   ------------------------------------------------------------



DefaultAssay(s4p_data)<-'SCT'

data.frame(
  orig.ident=s4p_data$orig.ident,
  cluster=s4p_data$fcluster,
  rpob=s4p_data@assays$Spatial@counts[c('rpoB'),],
  MAB0545=s4p_data@assays$Spatial@counts[c('MAB-0545'),]) %>% 
  filter(rpob+MAB0545>0) %>% 
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
  distinct() %>% 
  
  
  
  ggplot()+
  geom_violin(aes(x=gene,y=INOS,fill=gene),shape=21)+
  geom_boxplot(aes(x=gene,y=INOS),fill='white',width=0.4,shape=21)+
  facet_wrap(~orig.ident)+
  theme_classic()+
  ylab('INOS pathway levels')+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank())+
  scale_fill_manual(values=c('royalblue','#D95F02'),name='',labels=c('MABS+ lsr2-','MABS+ lsr2+'))

# Save Objects -------------------------------------------------------------



saveRDS(s4data,'Path/To/saved/object.rds')
saveRDS(s4p_data,'Path/To/saved/subsetted/object.rds')
