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


Load10X_Spatial('Path/to/cellranger/output/InfectedC1',filename = "filtered_feature_bc_matrix.h5",slice='SLICE')-> c_data

gene.cell.count = rowSums(GetAssayData(c_data,assay='Spatial')>0)
gene.above.threshold=names(gene.cell.count[gene.cell.count>2])
c_data=c_data[gene.above.theshold,]


Load10X_Spatial('Path/to/cellranger/output/InfectedD1',filename = "filtered_feature_bc_matrix.h5",slice='SLICE')-> d_data

gene.cell.count = rowSums(GetAssayData(d_data,assay='Spatial')>0)
gene.above.threshold=names(gene.cell.count[gene.cell.count>2])
d_data=d_data[gene.above.theshold,]


Load10X_Spatial('Path/to/cellranger/output/CTRLB1',filename = "filtered_feature_bc_matrix.h5",slice='SLICE')-> b_data

gene.cell.count = rowSums(GetAssayData(b_data,assay='Spatial')>0)
gene.above.threshold=names(gene.cell.count[gene.cell.count>2])
b_data=b_data[gene.above.theshold,]


Load10X_Spatial('Path/to/cellranger/output/CTRLA1',filename = "filtered_feature_bc_matrix.h5",slice='SLICE')-> a_data

gene.cell.count = rowSums(GetAssayData(a_data,assay='Spatial')>0)
gene.above.threshold=names(gene.cell.count[gene.cell.count>2])
a_data=a_data[gene.above.theshold,]


s4_data<-list(A1=a_data,B1=b_data,C1=c_data,D1=d_data)

# SCT normalization -------------------------------------------------------


s4_data$A1<-SCTransform(s4_data$A1,assay='Spatial',vst.flavor='v2')
s4_data$B1<-SCTransform(s4_data$B1,assay='Spatial',vst.flavor='v2')
s4_data$C1<-SCTransform(s4_data$C1,assay='Spatial',vst.flavor='v2')
s4_data$D1<-SCTransform(s4_data$D1,assay='Spatial',vst.flavor='v2')


# Perform deconvolution through RCTD --------------------------------------

s=c(A1='A1',B1='B1',D1='D1')

for (i in s){
  
  print(i)
  puck<-spacexr::SpatialRNA(s4_data[[i]]@images[[i]]@coordinates[,c('imagerow','imagecol')],s4_data[[i]]@assays$SCT@counts)
  
  
  spacexr::create.RCTD(puck,reference,max_cores = 2,test_mode = F, CELL_MIN_INSTANCE = 6)->myRCTD
  
  
  myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = 'full')
  
  
  
  results <- myRCTD@results
  
  norm_weights = as.matrix(sweep(results$weights, 1, rowSums(results$weights), '/'))
  norm_weights <- t(norm_weights)
  ncol <- length(colnames(s4_data[[i]])[!(colnames(s4_data[[i]]) %in% colnames(norm_weights))])
  if (ncol > 0) {
    tmp <- matrix(0, nrow = nrow(norm_weights), ncol = ncol)
    colnames(tmp) <- colnames(s4_data[[i]])[!(colnames(s4_data[[i]]) %in% colnames(norm_weights))]
    norm_weights <- cbind(norm_weights, tmp)[,colnames(s4_data[[i]])]
  }
  
  s4_data[[i]][['RCTD']]<-CreateAssayObject(data = norm_weights)
  
  
}

# Integrate slices ------------------------------------------------------------

features<-SelectIntegrationFeatures(object.list = s4_data)
s4_data<-PrepSCTIntegration(s4_data,anchor.features = features)
anchors<- FindIntegrationAnchors(object.list=s4_data,anchor.features=features,normalization.method = 'SCT')

s4data<-IntegrateData(anchorset = anchors)
s4data <- ScaleData(s4data, verbose = FALSE)
preprocessData(s4data,by='integrated')->s4data


DefaultAssay(s4data)<-'integrated'
s4data<-FindClusters(s4data,resolution = 0.22)
s4data$seurat_clusters<-s4data$integrated_snn_res.0.22


# Save Object -------------------------------------------------------------



saveRDS(s4data,'Path/To/saved/object.rds')
