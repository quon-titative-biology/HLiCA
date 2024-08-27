### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(colorspace)
library(cowplot)
library(SeuratDisk)
library(harmony)


source("scripts/00_pretty_plots.R")


obj <- readRDS(here("/media/redgar/Seagate Portable Drive/HLiCA/beta_annotation_objects/Mesenchyme_clean.rds"))
seu <- obj$seu
rm(obj)
gc()

DimPlot(seu, label=T)
DimPlot(seu, group.by="Beta.Annotation.SubLineage")



FeaturePlot(seu, "FOS")
FeaturePlot(seu, "JUN")

FeaturePlot(seu, "percent.mt")

DimPlot(seu, group.by = "STUDY")
DimPlot(seu, group.by = "suspension_type")

DotPlot(seu,features = c("FOS","JUN","EPCAM","KRT18","KRT7","KRT19","MUC5B","CXCL8","MSLN","LAMC2", "APOC3","APOA2","APOA1") )

de<-FindMarkers(seu, ident.1 = "5")



## Doublets

doublets<-seu@meta.data[which(seu$Beta.Annotation.SubLineage=="TBC"),]
doublets$Gamma.Annotation<-"doublets"

save(doublets, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/doublets_mesenchyme_meta.RData")


###############
## relabel
###############
cells_to_keep <- colnames(seu)[which(!(colnames(seu) %in% c(rownames(doublets))))]
seu_cleaned<-subset(seu, cells = cells_to_keep)
rm(seu)
gc()


## parameters match other runs by Jordan
seu_cleaned <- seu_cleaned %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% ScaleData()
seu_cleaned <- RunPCA(seu_cleaned, npcs = 60, verbose = FALSE)

covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')

seu_cleaned <- RunHarmony(seu_cleaned, group.by.vars = covariates, reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
seu_cleaned <- RunUMAP(seu_cleaned, reduction = "harmony", assay = "RNA", dims = 1:60)
seu_cleaned

seu_cleaned <- FindNeighbors(seu_cleaned, reduction = "harmony", dims = 1:60)
seu_cleaned <- FindClusters(seu_cleaned, resolution = 0.2)

plot_grid(DimPlot(seu_cleaned),DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage",label=T),DimPlot(seu_cleaned, group.by = "RNA_snn_res.0.1",label=T))


FeaturePlot(seu_cleaned, features = c("RELN") )
plot_grid(DimPlot(seu_cleaned, label=T),DotPlot(seu_cleaned,features = c("RELN","HHIP","COLEC10","VIPR1","LAMC3","SPON1","C7","FBLN1","RGS5","RGS6","MYH11") ))
plot_grid(DimPlot(seu_cleaned, label=T),DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage",label=T),DotPlot(seu_cleaned,features = c("RELN","HHIP","COLEC10","VIPR1","LAMC3","SPON1","C7","FBLN1","RGS5","RGS6","MYH11") ))

DimPlot(seu_cleaned, group.by="suspension_type",label=T)



unique(seu_cleaned$Beta.Annotation.SubLineage)
seu_cleaned$Gamma.Annotation<-as.character(seu_cleaned$seurat_clusters)

seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("1","3"))]<-"Hepatic Stellate Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("0","2","5","6","7","10","9"))]<-"Vascular Smooth Muscle Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("4"))]<-"Portal Fibroblast"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("8"))]<-"CUX2+ Hepatic Stellate Cell"




DimPlot(seu_cleaned, group.by="suspension_type")
head(de)

de<-FindMarkers(seu_cleaned, ident.1 = "8")
FeaturePlot(seu_cleaned, features = c("CUX2","SHANK2","LRP1B") )


de<-FindMarkers(seu_cleaned, ident.1 = "3")
de[which(de$avg_log2FC>0),][1:10,]
FeaturePlot(seu_cleaned, features = c("HP","SERPINA1","APOC3","ALB") )
plot_grid(DimPlot(seu_cleaned, label=T),DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage",label=T),DotPlot(seu_cleaned,features = c("RELN","HHIP","COLEC10","VIPR1","LAMC3","SPON1","C7","FBLN1","RGS5","RGS6","MYH11","CUX2","HP","SERPINA1","APOC3","ALB") ))


plot_grid(DimPlot(seu_cleaned),DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage",label=T),DimPlot(seu_cleaned, group.by = "Gamma.Annotation",label=T), DimPlot(seu_cleaned, group.by = "RNA_snn_res.0.2",label=T))


seu_cleaned$lineage_clusters<-seu_cleaned$RNA_snn_res.0.2


## save mesenchyme lineage object
save(seu_cleaned, file=here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_mesenchyme.RData"))

## save meta
mesenchyme<-seu_cleaned@meta.data
save(mesenchyme, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/mesenchyme_meta.RData")

