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


# Import data ----------------------------------------------------


Load10X_Spatial('Path/to/cellranger/output/Infected',filename = "filtered_feature_bc_matrix.h5",slice='SLICE')-> c_data


Load10X_Spatial('Path/to/cellranger/output/CTRL',filename = "filtered_feature_bc_matrix.h5",slice='SLICE')-> b_data
# SCT normalization -------------------------------------------------------

c_data<-SCTransform(c_data,vst.flavor = "v2",assay='Spatial')

b_data<-SCTransform(b_data,vst.flavor = "v2",assay='Spatial')


# Define Interested Areas ------------------------------------------------------------

read_delim('Path/to/spaceranger/barcodes/list',delim='\t',col_names = c('barcode','col','row'))->barcodes



## Infected

read_delim('path/to/list/of/selected/spots/Infected',delim=' ',col_names = c('col','row')) %>% 
  left_join(barcodes) %>% pull(barcode)->selected

# Example for parenchima
c_data@meta.data %>% rownames_to_column('spotID') %>% 
  mutate(Parenchima=str_remove(spotID,'-1')%in%selected) %>% 
  column_to_rownames('spotID')->c_data@meta.data


# Keep only the parenchima

Idents(c_data)<-'orig.ident'
cp_data<-subset(c_data,Parenchima==1)


## CTRL

read_delim('path/to/list/of/selected/spots/CTRL',delim=' ',col_names = c('col','row')) %>% 
  left_join(barcodes) %>% pull(barcode)->selected

# Example for parenchima
b_data@meta.data %>% rownames_to_column('spotID') %>% 
  mutate(Parenchima=str_remove(spotID,'-1')%in%selected) %>% 
  column_to_rownames('spotID')->b_data@meta.data


# Keep only the parenchima

Idents(b_data)<-'orig.ident'
bp_data<-subset(b_data,Parenchima==1)


# Perform deconvolution through RCTD --------------------------------------

## Infected
DefaultAssay(cp_data)<-'SCT'


puck<-spacexr::SpatialRNA(cp_data@images$C1@coordinates[,c('imagerow','imagecol')],cp_data@assays$SCT@counts)


## Use your rna_reference dataset from scRNAseq with cells definition


reference<-spacexr::Reference(rna_reference@assays$SCT@counts,rna_reference$cell_type)

spacexr::create.RCTD(puck,reference,max_cores = 2,test_mode = F, CELL_MIN_INSTANCE = 6)->myRCTD


myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = 'full')


results <- myRCTD@results

norm_weights = as.matrix(sweep(results$weights, 1, rowSums(results$weights), '/'))

norm_weights <- t(norm_weights)
ncol <- length(colnames(data)[!(colnames(data) %in% colnames(norm_weights))])
if (ncol > 0) {
  tmp <- matrix(0, nrow = nrow(norm_weights), ncol = ncol)
  colnames(tmp) <- colnames(data)[!(colnames(data) %in% colnames(norm_weights))]
  norm_weights <- cbind(norm_weights, tmp)[,colnames(data)]
}

cp_data[['RCTD']]<-CreateAssayObject(data = norm_weights)

DefaultAssay(cp_data)<-'RCTD'

## CTRL

DefaultAssay(bp_data)<-'SCT'


puck<-spacexr::SpatialRNA(bp_data@images$B1@coordinates[,c('imagerow','imagecol')],bp_data@assays$SCT@counts)


## Use your rna_reference dataset from scRNAseq with cells definition


reference<-spacexr::Reference(rna_reference@assays$SCT@counts,rna_reference$cell_type)

spacexr::create.RCTD(puck,reference,max_cores = 2,test_mode = F, CELL_MIN_INSTANCE = 6)->myRCTD


myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = 'full')


results <- myRCTD@results

norm_weights = as.matrix(sweep(results$weights, 1, rowSums(results$weights), '/'))

norm_weights <- t(norm_weights)
ncol <- length(colnames(data)[!(colnames(data) %in% colnames(norm_weights))])
if (ncol > 0) {
  tmp <- matrix(0, nrow = nrow(norm_weights), ncol = ncol)
  colnames(tmp) <- colnames(data)[!(colnames(data) %in% colnames(norm_weights))]
  norm_weights <- cbind(norm_weights, tmp)[,colnames(data)]
}

bp_data[['RCTD']]<-CreateAssayObject(data = norm_weights)

DefaultAssay(bp_data)<-'RCTD'



# Merge slices ------------------------------------------------------------



features<-SelectIntegrationFeatures(object.list = list(C1=cp_data,B1=bp_data))
anchors<- FindIntegrationAnchors(object.list=list(C1=cp_data,B1=bp_data),anchor.features=features)

scdata<-IntegrateData(anchorset = anchors)

scdata <- ScaleData(scdata, verbose = FALSE)

preprocessData(slide_merge)->slide_merge





SpatialDimPlot(slide_merge,'SCT_snn_res.0.3',alpha=0.8,crop=F,pt.size.factor =1)


slide_merge$seurat_clusters<-slide_merge$SCT_snn_res.0.4



# Umap and plots ----------------------------------------------------------



SpatialDimPlot(scdata,'seurat_clusters',alpha=0.8,crop=F,pt.size.factor =1,combine=F)->p



wrap_plots(p[[2]],p[[1]])


scdata$Infection=ifelse(str_detect(names(scdata$orig.ident),'_1'),'Chronic Infected','CTRL')

wrap_plots(
  DimPlot(scdata,group.by = 'seurat_clusters')+
    ggtitle('Clusters')+#DarkTheme()+
    scale_color_manual(name='',values = gg_color_hue(11)),
  DimPlot(scdata,group.by = c('Infection'))+
    scale_color_manual(values=c('red2','royalblue')))



# Infected Slice analysis -------------------------------------------------

## Get mabs spot

DefaultAssay(cp_data)<-'Spatial'

sum_counts <- cp_data@assays$Spatial@counts['rpoB', ] + cp_data@assays$Spatial@counts['MAB-0545', ]

cp_data$mabs<-sum_counts


cp_data$mabs_spot<-ifelse(cp_data$mabs==0,0,1)



## Get Distance from MABS spot

mtx0 = cp_data@images[['C1']]@coordinates
mtx0 <- cbind(mtx0, t(cp_data@assays$RCTD@data))
mtx0 <- cbind(mtx0, cp_data@meta.data)
mtx0$Mgp <- cp_data@assays$SCT@data['Mgp',]
airCells <-  data.matrix(mtx0[mtx0['mabs_spot'] == 1, c("imagerow", "imagecol")])
cell <- data.matrix(mtx0[c("imagerow", "imagecol", rownames(cp_data@assays$RCTD@data))])
#cell.dist <- computeDist(cell, airCells)
nAir = nrow(airCells)
helper = matrix(rep(1,nrow(cell)), ncol=1)

outMtx = matrix(nrow = nrow(subCell), ncol = nAir)
rownames(outMtx) = rownames(subCell)
colnames(outMtx) = rownames(airCells)

for (i in 1:nAir) {
  temp0 = subCell[,c(1,2)] - helper %*% airCells[i,]
  outMtx[,i] = sqrt(temp0[,1]^2 + temp0[,2]^2)
}
min.dis = apply(outMtx,1,min)
return(min.dis)


cell %>%as.data.frame() %>%  rownames_to_column('spotID') %>% mutate(a='a') %>% 
  left_join(airCells %>% as.data.frame() %>%  rownames_to_column('spotIDMabs') %>% 
              transmute(spotIDMabs=spotIDMabs,row=imagerow,col=imagecol,a='a')) %>% 
  mutate(dist=sqrt((row-imagerow)^2+(col-imagecol)^2)) %>% 
  gather(Type,perc,`Alveolar epithelial cells`:Erythrocytes)->mabdist


## Get Near mabs spot


mabdist %>% 
  select(-Type) %>% distinct() %>% 
  filter(dist<200) %>% 
  select(spotID) %>% distinct() %>% pull(spotID)->neraMabs



cp_data$near_mabs<-rownames(cp_data@meta.data)%in%neraMabs



## Get Only spots Mabs postives and neighbours

c1p_subs<-subset(cp_data,near_mabs==1)


preprocessData(c1p_subs)->c1p_subs
SpatialDimPlot(c1p_subs,group.by = paste('SCT_snn_res.',seq(0.1,1,0.1),sep=''),combine = F)

c1p_subs$seurat_clusters<-c1p_subs$SCT_snn_res.0.4



Idents(c1p_subs)<-'seurat_clusters'

## Map of clusters near Mabs

c1p_subs$FClusters<-factor(c1p_subs$Clusters,levels=str_wrap(c('Upper Airways','Upper Airways - Submucosal','Lower Airways - Parenchima','Lower Airways - Granuloma-like structure')))

SpatialDimPlot(c1p_subs,crop=F,pt.size.factor = 1,image.alpha = 0.7,group.by = 'FClusters')+
  theme(legend.title = element_blank())

## Barplots

data.frame(cluster=c1p_subs$SCT_snn_res.0.4,
           rpob=c1p_subs@assays$Spatial@counts[c('rpoB'),],
           MAB0545=c1p_subs@assays$Spatial@counts[c('MAB-0545'),]) %>% 
  filter(rpob+MAB0545>0) %>% 
  gather(gene,count,-cluster) %>% group_by(gene,cluster) %>% 
  summarise(c=sum(count)) %>% 
  ungroup() %>% 
  mutate(cluster=factor(cluster,levels=c(2,3,0,1))) %>% 
  mutate(gene=factor(gene,levels=c('rpob','MAB0545'))) %>% 
  ggplot(aes(y=c,x=cluster))+
  geom_bar(stat='identity',position='dodge',col='black',aes(fill=gene))+
  theme_classic()+
  scale_fill_manual(values=rev(c('#D95F02','#1B9E77')),name='',labels=rev(c('lsr2','rpoB')))+
  scale_x_discrete(labels=str_wrap(c('Upper Airways','Upper Airways - Submucosal','Lower Airways - Parenchima','Lower Airways - Granuloma-like structure'),15))+
  theme(axis.text.x = element_text(angle=45,vjust=0.6),
        axis.title.x = element_blank())+
  ylab('UMI counts')+
  geom_text(data=chisq,aes(y=45,x=0.5,label=paste('Chi-squared test:',signif(p,1))),hjust=0)


## Percentage Barplots + Inos


nos2markers<-c('Acp5','Cav1','Cyb5b','Cyb5r3','Cyp1b1','Ddah1','H2-M3','Hsp90aa1','Ifng','Klrk1','Mtarc1','Mtarc2','Nos1','Nos1','Nos1','Nos1','Nos2','Nos2','Nos2','Nos2','Nos2','Nos3','Nos3','Nos3','Nos3','P2rx4','P2rx4','Rora','Slc7a2','Slc7a6','Spr','Ticam1','Tlr2','Tlr2','Tlr4')


AddModuleScore(c1p_subs,list(c(nos2markers)),name='INOS')->c1p_subs


data.frame(cluster=c1p_subs$SCT_snn_res.0.4,
           rpob=c1p_subs@assays$Spatial@counts[c('rpoB'),],
           MAB0545=c1p_subs@assays$Spatial@counts[c('MAB-0545'),]) %>% 
  filter(rpob+MAB0545>0) %>% group_by(cluster) %>% add_count() %>% 
  gather(gene,count,-cluster,-n) %>% group_by(cluster,gene) %>% 
  summarise(c=sum(count)) %>% mutate(cc=sum(c)) %>%distinct() %>% 
  left_join(data.frame(cluster=c1p_subs$SCT_snn_res.0.4,
                       INOS=c1p_subs$INOS1) %>% 
              group_by(cluster) %>% summarise(Inos=median(INOS))) %>% 
  ungroup() %>% 
  mutate(cluster=factor(cluster,levels=c(2,3,0,1))) %>% 
  ggplot()+
  geom_col(aes(x=cluster,y=c,fill=gene),position = 'fill',col='black',width=0.8) +
  geom_line(aes(x=as.integer(cluster),y=(Inos+0.3693719)/(0.4579331+0.3693719)),data=~select(.,-gene) %>% distinct(),lty=2)+
  # geom_violin(aes(x=cluster,y=(INOS+0.3693719)/(0.4579331+0.3693719)),data=data.frame(cluster=c1p_subs$SCT_snn_res.0.4,
  #                    INOS=c1p_subs$INOS1),width=0.3)+
  geom_boxplot(aes(x=cluster,y=(INOS+0.3693719)/(0.4579331+0.3693719)),data=data.frame(cluster=c1p_subs$SCT_snn_res.0.4,
                                                                                       INOS=c1p_subs$INOS1),width=0.2,lwd=1)+
  scale_fill_manual(values=c('#D95F02','#1B9E77'),name='',labels=c('lsr2','rpoB'))+
  scale_y_continuous(sec.axis =sec_axis(~(.x*(0.4579331+0.3693719))-0.3693719,name=str_wrap(str_to_title('nitric oxide biosynthetic process pathway levels'),30)),name='UMI counts distribution')+
  theme_classic()+
  scale_x_discrete(labels=str_wrap(c('Upper Airways','Upper Airways - Submucosal','Lower Airways - Parenchima','Lower Airways - Granuloma-like structure'),15))+
  theme(axis.text.x = element_text(angle=45,vjust=0.6),
        axis.title.x = element_blank())


## Map of INOS

SpatialFeaturePlot(c1p_subs,'INOS1',crop=F,pt.size.factor = 1,image.alpha = 0.7) +
  scale_fill_gradient(low = 'lightgreen',high='red',name=str_wrap(str_to_title('nitric oxide biosynthetic process'),20))


## Distances graphs

dist_matrix<- mabdist %>% 
  filter(dist<1000) %>% 
  mutate(dist=(round(dist,0))) %>% 
  select(spotIDMabs,dist,perc,Type) %>% 
  group_by(spotIDMabs,dist,Type) %>% 
  summarise(perc=mean(perc))


lapply(names(table(c1p_subs$Clusters)),function(x){
  set.seed(0)
  dist_matrix %>% 
    left_join(c1p_subs@meta.data %>% rownames_to_column('spotIDMabs') %>% select(spotIDMabs,cluster=Clusters)) %>% 
    #filter(Type%in%colnames(cell)) %>% 
    filter(Type%in%c('B cells','Dendritic cells','Macrophages','Monocytes','Neutrophils','T cells')) %>% 

    filter(cluster==x) %>% 
    ggplot(aes(dist, perc, color=Type)) + 
    geom_smooth(method = "loess") + scale_color_discrete(type=sample(palette)) + 
    facet_wrap(~cluster)+theme_bw()+

    theme(strip.background=element_rect(colour="black",
                                        fill="white"),
          strip.text = element_text(size=15),
          legend.title = element_blank(),
          axis.title = element_text(size=15))+
    scale_y_continuous(labels = scales::percent)+
    ylab('Proportion of cell types')+xlab(expression(paste('Distance from M.abs probes(',mu,'m)')))
  
})

