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


obj <- readRDS(here("data/Myeloid_Cell_clean.rds"))
seu <- obj$seu
rm(obj)
gc()

DimPlot(seu)
DimPlot(seu, group.by="Beta.Annotation.SubLineage")



## pDC
pDC<-seu@meta.data[which(seu$Beta.Annotation.SubLineage=="pDCs"),]
pDC$Gamma.Annotation<-"pDC"
save(pDC, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/pDC_meta.RData")

## Doublets
## LSEC doublets are with Kupffer cells
d10x.combined_KC<-subset(seu, subset = Beta.Annotation.SubLineage %in% c("Kupffer Cells"))
d10x.combined_KC <- RunPCA(d10x.combined_KC, npcs = 30, verbose = FALSE)
d10x.combined_KC <- RunUMAP(d10x.combined_KC, reduction = "pca", dims = 1:30)
d10x.combined_KC <- FindNeighbors(d10x.combined_KC, reduction = "pca", dims = 1:30)
d10x.combined_KC <- FindClusters(d10x.combined_KC, resolution = 0.2)

DimPlot(d10x.combined_KC, label=T)
DimPlot(d10x.combined_KC, group.by = "Beta.Annotation.SubLineage")

FeaturePlot(d10x.combined_KC, "CALCRL")
FeaturePlot(d10x.combined_KC, "RAMP2")

LSEC_doublets<-d10x.combined_KC@meta.data[which(d10x.combined_KC$seurat_clusters=="11"),]
tcell_doublets<-seu@meta.data[which(seu$Beta.Annotation.SubLineage=="T Cells?"),]

doublets<-rbind(tcell_doublets, LSEC_doublets)
doublets$Gamma.Annotation<-"doublets"

save(doublets, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/doublets_myeloid_meta.RData")


  
### remove and recluster
load("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/doublets_myeloid_meta.RData")
load("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/pDC_meta.RData")

cells_to_keep <- colnames(seu)[which(!(colnames(seu) %in% c(rownames(doublets), rownames(pDC))))]
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



### review markers discussed
DimPlot(seu_cleaned, label=T)
plot_grid(DimPlot(seu_cleaned, label=T), DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage"), rel_widths = c(1,1.25))
save_plts(plot_grid(DimPlot(seu_cleaned, label=T), DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage"), rel_widths = c(1,1.25)), "recluster_myeloids", h=5,w=15)

DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage")+scale_color_manual(values = c(rep("grey", 2), "red",rep("grey", 5),"blue"))
DimPlot(seu_cleaned, group.by = "STUDY")

FeaturePlot(seu_cleaned, features = "MARCO")
FeaturePlot(seu_cleaned, features = "TREM2")
FeaturePlot(seu_cleaned, features = "LAMP3")
FeaturePlot(seu_cleaned, features = c("BAG3"))
FeaturePlot(seu_cleaned, features = c("HLA-DRA"))
FeaturePlot(seu_cleaned, features = c("CLEC9A", "CLEC10A"))
FeaturePlot(seu_cleaned, features = c("MAMLD1", "IGSF21"))

FeaturePlot(seu_cleaned, features = c("TREM2","FOLR2","EREG","VCAN"))

FeaturePlot(seu_cleaned, features = c("MARCO", "C1QC", "CD5L", "LYZ", "CCR2","MRC1"))

DotPlot(seu_cleaned,features = c("TREM2","FOLR2","EREG","VCAN","LAMP3","PLEK2", "CLEC10A", "FCER1A", "IGSF21","CTLA4","MARCO", "C1QC", "CD5L") )
DotPlot(seu_cleaned,features = c("CCR7", "LAMP3", "BIRC3", "CD274","PLEK2", "CLEC10A", "FCER1A","CTLA4") )
DotPlot(seu_cleaned, features = c("MARCO", "C1QC", "CD5L", "LYZ", "CCR2", "MRC1"))

# ## label neutrophils
FeaturePlot(seu_cleaned, features=c("FCGR3B","CSF3R","S100P","MME"))




de<-FindMarkers(seu_cleaned, ident.1 = "0", ident.2="2")
de<-FindMarkers(seu_cleaned, ident.1 = "0", ident.2="1")
de<-FindMarkers(seu_cleaned, ident.1 = "17", ident.2="1")

VlnPlot(seu_cleaned, features = c("FCGR3B", "RMC1", "nFeature_RNA"), pt.size = 0)

DotPlot(seu_cleaned, features = c("MARCO", "C1QC", "CD5L", "LYZ", "FCGR3B","CMTM2","HBB"))


###############
## relabel
###############
unique(seu_cleaned$Beta.Annotation.SubLineage)
seu_cleaned$Gamma.Annotation<-as.character(seu_cleaned$seurat_clusters)

seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("11"))]<-"TREM2+ Macrophages"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("8"))]<-"Neutrophils"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("18"))]<-"migDCs"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("14"))]<-"DCmac"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("16"))]<-"BAG3+ Monocytes"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("9","15"))]<-"MAMLD1+ Trans Monocytes"

## previous label correct
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("13"))]<-"Cycling Cells"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("10","19"))]<-"Type 1 cDCs"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("5","6"))]<-"Type 2 cDCs"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("0","7"))]<-"Kupffer Cells"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("4"))]<-"Non-classical Monocytes"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("1","2","3","17"))]<-"Classical Monocytes"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("12"))]<-"Activated Monocytes"

DimPlot(seu_cleaned, group.by = "Gamma.Annotation",label=T)
DimPlot(seu_cleaned, group.by = "RNA_snn_res.0.6",label=T)

seu_cleaned$lineage_clusters<-seu_cleaned$RNA_snn_res.0.6

## save myeloid lineage object
save(seu_cleaned, file=here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_cleaned_myeloid.RData"))

## save meta
myeloid<-seu_cleaned@meta.data
save(myeloid, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/myeloid_meta.RData")


# ## as h5ad
# load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_cleaned_myeloid.RData"))
# DefaultAssay(seu_cleaned)<-"RNA"
# SaveH5Seurat(seu_cleaned, filename = here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_cleaned_myeloid.h5Seurat"), overwrite=T)
# Convert(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_cleaned_myeloid.h5Seurat"), dest = "h5ad",overwrite=T)

DotPlot(seu_cleaned, group.by = "RNA_snn_res.0.6",features = c("MARCO", "C1QC", "CD5L", "LYZ", "CCR2", "MRC1","CD1C","CLEC10A","FCER1A","HLA-DQA1"))
