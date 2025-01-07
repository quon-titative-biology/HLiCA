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
source("scripts/00_fanciest_UMAP.R")
source("scripts/colour_palette.R")



obj <- readRDS(here("/media/redgar/Seagate Portable Drive/HLiCA/beta_annotation_objects/Lymphocyte_clean.rds"))
seu <- obj$seu
rm(obj)
gc()

DimPlot(seu)
DimPlot(seu, group.by="Beta.Annotation.SubLineage")



## Doublets

doublets<-seu@meta.data[which(seu$Beta.Annotation.SubLineage%in%c("Stressed Cells","Doublet","Poor Quality")),]
doublets$Gamma.Annotation<-"doublets"

save(doublets, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/doublets_lymphocytes_meta.RData")


######################### ON COMPUTE CAN  
obj <- readRDS("/home/redgar25/scratch/Lymphocyte_clean.rds")
seu <- obj$seu
rm(obj)
gc()


## Doublets

doublets<-seu@meta.data[which(seu$Beta.Annotation.SubLineage%in%c("Stressed Cells","Doublet","Poor Quality")),]
doublets$Gamma.Annotation<-"doublets"


cells_to_keep <- colnames(seu)[which(!(colnames(seu) %in% c(rownames(doublets))))]
seu_cleaned<-subset(seu, cells = cells_to_keep)
rm(seu)
gc()


## parameters match other runs by Jordan
seu_cleaned <- seu_cleaned %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% ScaleData()
seu_cleaned <- RunPCA(seu_cleaned, npcs = 40, verbose = FALSE)

covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')

seu_cleaned <- RunHarmony(seu_cleaned, group.by.vars = covariates, reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
seu_cleaned <- RunUMAP(seu_cleaned, reduction = "harmony", assay = "RNA", dims = 1:40)
seu_cleaned

seu_cleaned <- FindNeighbors(seu_cleaned, reduction = "harmony", dims = 1:40)
seu_cleaned <- FindClusters(seu_cleaned, resolution = 0.3)
seu_cleaned <- FindClusters(seu_cleaned, resolution = 0.4)
seu_cleaned <- FindClusters(seu_cleaned, resolution = 0.6)


save(seu_cleaned, file="/home/redgar25/scratch/tmp_lymphocytes_seurat.RData")
######################### ON COMPUTE CAN  





load("/media/redgar/Seagate Portable Drive/HLiCA/beta_annotation_objects/tmp_lymphocytes_seurat.RData")

### review markers discussed
DimPlot(seu_cleaned, label=T)
plot_grid(DimPlot(seu_cleaned, group.by="RNA_snn_res.0.4",label=T), DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage"), rel_widths = c(1,1.25))
save_plts(plot_grid(DimPlot(seu_cleaned,group.by="RNA_snn_res.0.4", label=T), DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage"), rel_widths = c(1,1.25)), "recluster_lymphocytes", h=5,w=15)

table(seu_cleaned$RNA_snn_res.0.4, seu_cleaned$Beta.Annotation.SubLineage)

DimPlot(seu_cleaned, group.by = "RNA_snn_res.0.6", label=T)
DimPlot(seu_cleaned, group.by = "RNA_snn_res.0.4", label=T)
DimPlot(seu_cleaned, group.by = "RNA_snn_res.0.3", label=T)


FeaturePlot(seu_cleaned, features = c("PDCD1","HAVCR2","CTLA4"))
DotPlot(seu_cleaned, group.by = "RNA_snn_res.0.4", features = c("PDCD1","HAVCR2","CTLA4"))


###############
## relabel
###############
unique(seu_cleaned$Beta.Annotation.SubLineage)
seu_cleaned$Gamma.Annotation<-as.character(seu_cleaned$RNA_snn_res.0.4)

seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("0","10","11","13"))]<-"CD8 T Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("5","12"))]<-"Exhausted CD8 T Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("1"))]<-"Bright NK Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("2"))]<-"MAIT T Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("3"))]<-"Dim NK Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("4"))]<-"Helper T Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("6"))]<-"Plasma B Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("7"))]<-"B Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("8"))]<-"Cycling"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("9"))]<-"Regulatory T Cell"


plasma<-subset(seu_cleaned, subset = Gamma.Annotation == 'Plasma B Cell')
plasma <- plasma %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% ScaleData()
plasma <- RunPCA(plasma, npcs = 40, verbose = FALSE)
plasma <- RunUMAP(plasma, reduction = "harmony", assay = "RNA", dims = 1:40)
plasma <- FindNeighbors(plasma, reduction = "harmony", dims = 1:40)
plasma <- FindClusters(plasma, resolution = 0.2)

DimPlot(plasma)
FeaturePlot(plasma, features = c("IGHG1","IGHA1"))
DotPlot(plasma, features = c("IGHG1","IGHA1"))

plasma$Gamma.Annotation[which(plasma$seurat_clusters%in%c("1"))]<-"IgA B cells"
plasma$Gamma.Annotation[which(plasma$seurat_clusters%in%c("0"))]<-"IgG B cells"
plasma$Gamma.Annotation[which(plasma$seurat_clusters%in%c("2","3"))]<-"Plasma B Cell"

DimPlot(plasma, group.by = "Gamma.Annotation")
DotPlot(plasma, group.by = "Gamma.Annotation", features = c("IGHG1","IGHA1"))

IGG<-colnames(subset(plasma, subset = Gamma.Annotation == 'IgG B cells'))
IGA<-colnames(subset(plasma, subset = Gamma.Annotation == 'IgA B cells'))

seu_cleaned$Gamma.Annotation[which(colnames(seu_cleaned) %in% IGG )]<-"IgG B cells"
seu_cleaned$Gamma.Annotation[which(colnames(seu_cleaned) %in% IGA )]<-"IgA B cells"

DimPlot(seu_cleaned, group.by = "Gamma.Annotation")

Idents(seu_cleaned)<-"RNA_snn_res.0.4"
#de<-FindMarkers(seu_cleaned, ident.1 = "10")

# 
# VlnPlot(seu_cleaned, features = c("FCGR3B", "RMC1", "nFeature_RNA"), pt.size = 0)
# 
# DotPlot(seu_cleaned, features = c("MARCO", "C1QC", "CD5L", "LYZ", "FCGR3B","CMTM2","HBB"))


seu_cleaned$lineage_clusters<-seu_cleaned$RNA_snn_res.0.4

## save lineage object
save(seu_cleaned, file=here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_cleaned_lymphocytes.RData"))

## save meta
lymphocyte<-seu_cleaned@meta.data
save(lymphocyte, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/lymphocyte_meta.RData")




#### Fancy UMAP

load(("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_cleaned_lymphocytes.RData"))
UMAP = as.data.frame(Embeddings(seu_cleaned, reduction = "umap"))

meta<-seu_cleaned@meta.data
rm(seu_cleaned)
gc()

meta$Gamma.Annotation[which(meta$Gamma.Annotation=="Cycling Cells")]<-"Cycling"

fanciest_UMAP(meta, UMAP)
save_plts(fanciest_UMAP(meta, UMAP), "lymphocytes_fancy", w=5, h=3.5)

fanciest_UMAP(meta, UMAP, rnd_col = T)
save_plts(fanciest_UMAP(meta, UMAP,rnd_col = T), "lymphocytes_fancy_rndcol", w=5, h=3.5)
