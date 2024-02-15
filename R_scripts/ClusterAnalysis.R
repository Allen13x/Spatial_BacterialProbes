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




# Load Object -------------------------------------------------------------

s4data<-readRDS('Path/To/saved/object.rds')




# Label and order Transcriptional Clusters --------------------------------



s4data$Infection=ifelse(s4data$orig.ident%in%c('C1','D1'),'Chronic Infected','CTRL')
s4data@meta.data %>% 
  rownames_to_column('spot') %>% 
  group_by(seurat_clusters) %>%
  mutate(nn=n()) %>% 
  group_by(Infection,seurat_clusters) %>% 
  mutate(n=n(),
         c_id=case_when((Infection=='CTRL'& n>0.9*nn)~'C',
                        (Infection!='CTRL'& n>0.9*nn)~'I',
                        (Infection=='CTRL'& n<0.1*nn)~'I',
                        (Infection!='CTRL'& n<0.1*nn)~'C',
                        T~'S')) %>%
  ungroup() %>% 
  group_by(c_id) %>% 
  mutate(c_n=as.integer(factor(seurat_clusters)),
         c_id=paste0(c_id,c_n)) %>% 
  ungroup() %>% 
  mutate(Clusters=(factor(c_id,levels=c(paste0('S',c(1:5)),'C1','I1','I2')))) %>% 
  column_to_rownames('spot')->s4data@meta.data



# Cluster markers ---------------------------------------------------------

DefaultAssay(s4data)<-'SCT'
Idents(s4data)<-'Clusters'
s4data<-PrepSCTFindMarkers(s4data)
FindAllMarkers(s4data)->s4data_markers





# Umap and plots ----------------------------------------------------------

Idents(s4data)<-'Clusters'

Cluster_palette=c('#d67d0f',
                  '#60b038',
                  '#9862d1',
                  '#decc21',
                  '#62d1c6',
                  'grey',
                  '#e60b0b',
                  '#ed555f')

names(Cluster_palette)=s

SpatialDimPlot(s4data,'Clusters',alpha=1,crop=F,pt.size.factor =1.5,combine=F,image.alpha = 0,stroke=0)->p


wrap_plots(p[[4]]+
             scale_fill_manual(values=Cluster_palette[-6]),
           p[[3]]+
             scale_fill_manual(values=Cluster_palette[-6]),
           p[[1]]+
             scale_fill_manual(values=Cluster_palette[-6]),
           p[[2]]+
             scale_fill_manual(values=Cluster_palette))



s4data$Infection=ifelse(s4data$orig.ident%in%c('C1','D1'),'Chronic Infection','Ctrl')

wrap_plots(
  DimPlot(s4data,group.by = 'orig.ident')+
    ggtitle('Samples')+
    scale_color_manual(values=c('turquoise','royalblue','red2','orange'),name=''),
  DimPlot(s4data,group.by = 'Clusters')+
    ggtitle('Clusters')+
    scale_color_manual(name='',values = Cluster_palette),
  DimPlot(s4data,group.by = c('Infection'))+
    scale_color_manual(name='',values=c('red2','royalblue')))





# Heatmap -----------------------------------------------------------------

s4data_markers %>% 
  group_by(cluster) %>% 
  filter(!str_detect(gene,'^Hb')) %>% 
  select(gene,avg_log2FC,cluster) %>% 
  spread(cluster,avg_log2FC) %>% 
  filter(gene%in%(s4data_markers %>%
                    filter(!str_detect(gene,'^Hb')) %>% 
                    group_by(cluster) %>%
                    slice_max(avg_log2FC,n=7) %>% pull(gene))) %>% 
  replace(is.na(.),0) %>% 
  column_to_rownames('gene')->mm


pheatmap(mm[s4data_markers %>%
              filter(!str_detect(gene,'^Hb')) %>% 
              group_by(cluster) %>%
              slice_max(avg_log2FC,n=7) %>% pull(gene),],
         cluster_rows = F,cluster_cols = F,
         annotation_col = data.frame(Cluster=names(Cluster_palette),row.names = names(Cluster_palette)),
         annotation_colors = list(Cluster=Cluster_palette),
         breaks = seq(-1,1,length.out=101),
         gaps_col = seq(0, 8),
         labels_row=s4data_markers %>%
           filter(!str_detect(gene,'^Hb')) %>% 
           group_by(cluster) %>%
           slice_max(avg_log2FC,n=7) %>% pull(gene),
         show_colnames = F,
         scale='none',
         color = colorRampPalette((c('magenta','magenta','black','yellow','yellow')))(100))




# Cellular Composition Cluster--------------------------------------


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
  mutate(perc=filter_lims(perc)) %>% 
  ggplot(aes(x=var,y=perc)) +
  geom_boxplot(aes(fill=var),na.rm = T,coef=7,lwd=0.5,show.legend = T)+
  theme_classic()+
  ggtitle(str_to_title(str_replace_all(var,'_',' ')))+
  theme(legend.title = element_blank(),
        legend.position = 'bottom',
        axis.ticks.x = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=10),
        axis.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
facet_wrap(~splitvar,scales='fixed',ncol = 4)




# Save Object -------------------------------------------------------------



saveRDS(s4data,'Path/To/saved/object.rds')

