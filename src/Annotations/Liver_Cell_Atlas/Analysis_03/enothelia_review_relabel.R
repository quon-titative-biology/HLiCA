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


obj <- readRDS(here("/media/redgar/Seagate Portable Drive/HLiCA/beta_annotation_objects/Endothelia_clean.rds"))
seu_cleaned <- obj$seu
rm(obj)
gc()

DimPlot(seu_cleaned)
DimPlot(seu_cleaned, group.by="Beta.Annotation.SubLineage")








###############
## relabel
###############
unique(seu_cleaned$Beta.Annotation.SubLineage)
seu_cleaned$Gamma.Annotation<-as.character(seu_cleaned$seurat_clusters)

DimPlot(seu_cleaned, group.by="Gamma.Annotation")

## previous label correct
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("1","3","5","6","7"))]<-as.character(seu_cleaned$Beta.Annotation.SubLineage[which(seu_cleaned$Gamma.Annotation%in%c("1","3","5","6","7"))])
DimPlot(seu_cleaned, group.by = "Gamma.Annotation",label=T)

seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("0"))]<-"Periportal LSEC"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("2"))]<-"Central Venous LSEC"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("4"))]<-"Vascular Endothelial Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("6"))]<-"Lymphatic Endothelial Cell"


plot_grid(DimPlot(seu_cleaned, group.by = "Gamma.Annotation",label=T),DimPlot(seu_cleaned, group.by = "seurat_clusters",label=T))

seu_cleaned$lineage_clusters<-seu_cleaned$seurat_clusters

## save endo lineage object
save(seu_cleaned, file=here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_cleaned_endo.RData"))

## save meta
endo<-seu_cleaned@meta.data
save(endo, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/endo_meta.RData")



#### Fancy UMAP

load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_cleaned_endo.RData"))
UMAP = as.data.frame(Embeddings(seu_cleaned, reduction = "umap"))

meta<-seu_cleaned@meta.data
rm(seu_cleaned)
gc()

meta$Gamma.Annotation[which(meta$Gamma.Annotation=="Cycling Cells")]<-"Cycling"

colnames(UMAP)[which(colnames(UMAP)=="UMAP_1")]<-"umap_1"
colnames(UMAP)[which(colnames(UMAP)=="UMAP_2")]<-"umap_2"

fanciest_UMAP(meta, UMAP)
save_plts(fanciest_UMAP(meta, UMAP), "endothelia_fancy", w=5, h=3.5)

fanciest_UMAP(meta, UMAP, rnd_col = T)
save_plts(fanciest_UMAP(meta, UMAP,rnd_col = T), "endothelia_fancy_rndcol", w=5, h=3.5)




