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



obj <- readRDS(here("/media/redgar/Seagate Portable Drive/HLiCA/beta_annotation_objects/Cholangiocyte_clean.rds"))
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

save(doublets, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/doublets_cholangiocytes_meta.RData")


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
seu_cleaned <- FindClusters(seu_cleaned, resolution = 0.6)

plot_grid(DimPlot(seu_cleaned),DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage",label=T))

FeaturePlot(seu_cleaned, features = c("LAMC2") )
plot_grid(DimPlot(seu_cleaned, label=T),DotPlot(seu_cleaned,features = c("FOS","JUN","EPCAM","KRT18","KRT7","KRT19","MUC5B","CXCL8","MSLN","LAMC2", "APOC3","APOA2","APOA1") ))

unique(seu_cleaned$Beta.Annotation.SubLineage)
seu_cleaned$Gamma.Annotation<-as.character(seu_cleaned$seurat_clusters)

seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("0","1","3"))]<-"Small ApoLipo"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("6"))]<-"LAMC2+ Small"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("7"))]<-"CXCL8+ Small Keratin"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("2","8","9","4"))]<-"Small Keratin"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("5"))]<-"Large Mucus Secreting"


de<-FindMarkers(seu_cleaned, ident.1 = "4",ident.2 = "2")
FeaturePlot(seu_cleaned, features = c("RPS3A", "RPS27A","MUC5B") )
DimPlot(seu_cleaned, group.by="suspension_type")
de<-FindMarkers(seu_cleaned, ident.1 = "7")


plot_grid(DimPlot(seu_cleaned),DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage",label=T),DimPlot(seu_cleaned, group.by = "Gamma.Annotation",label=T), DimPlot(seu_cleaned, group.by = "RNA_snn_res.0.2",label=T))


seu_cleaned$lineage_clusters<-seu_cleaned$RNA_snn_res.0.6


## save cholangiocyte lineage object
save(seu_cleaned, file=here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_cholangiocyte.RData"))

## save meta
cholangiocyte<-seu_cleaned@meta.data
save(cholangiocyte, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/cholangiocyte_meta.RData")




#### Fancy UMAP

load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_cholangiocyte.RData"))
UMAP = as.data.frame(Embeddings(seu_cleaned, reduction = "umap"))

meta<-seu_cleaned@meta.data
rm(seu_cleaned)
gc()

meta$Gamma.Annotation[which(meta$Gamma.Annotation=="Cycling Cells")]<-"Cycling"

fanciest_UMAP(meta, UMAP)
save_plts(fanciest_UMAP(meta, UMAP), "cholangiocyte_fancy", w=5, h=3.5)

fanciest_UMAP(meta, UMAP, rnd_col = T)
save_plts(fanciest_UMAP(meta, UMAP,rnd_col = T), "cholangiocyte_fancy_rndcol", w=5, h=3.5)


