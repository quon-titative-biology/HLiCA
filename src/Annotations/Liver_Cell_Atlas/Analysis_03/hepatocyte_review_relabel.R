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
#library(SeuratDisk)
library(harmony)


source("scripts/00_pretty_plots.R")
source("scripts/00_fanciest_UMAP.R")
source("scripts/colour_palette.R")


######################### ON COMPUTE CAN  
obj <- readRDS("/home/redgar25/scratch/Hepatocyte_clean.rds")
seu <- obj$seu
rm(obj)
gc()


## Doublets

remove<-seu@meta.data[which(seu$Beta.Annotation.SubLineage%in%c("Mito+","Potentially Technical")),]
remove$Gamma.Annotation<-"removed"
save(remove, file="/home/redgar25/scratch/remove_hepatocytes_meta.RData")


cells_to_keep <- colnames(seu)[which(!(colnames(seu) %in% c(rownames(remove))))]
seu_cleaned<-subset(seu, cells = cells_to_keep)
rm(seu)
gc()


## parameters match other runs by Jordan
seu_cleaned <- seu_cleaned %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData()
seu_cleaned <- RunPCA(seu_cleaned, npcs = 20, verbose = FALSE)

covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')

seu_cleaned <- RunHarmony(seu_cleaned, group.by.vars = covariates, reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
seu_cleaned <- RunUMAP(seu_cleaned, reduction = "harmony", assay = "RNA", dims = 1:20)
seu_cleaned

seu_cleaned <- FindNeighbors(seu_cleaned, reduction = "harmony", dims = 1:20)
seu_cleaned <- FindClusters(seu_cleaned, resolution = 0.3)
seu_cleaned <- FindClusters(seu_cleaned, resolution = 0.4)
seu_cleaned <- FindClusters(seu_cleaned, resolution = 0.6)


save(seu_cleaned, file="/home/redgar25/scratch/tmp_hepatocytes_seurat.RData")


ggsave(plot_grid(DimPlot(seu_cleaned,group.by="RNA_snn_res.0.3", label=T), DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage"), rel_widths = c(1,1.25)), file="/home/redgar25/scratch/recluster_hepatocytes.jpeg", h=5,w=15)
ggsave(plot_grid(DimPlot(seu_cleaned,group.by="RNA_snn_res.0.4", label=T), DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage"), rel_widths = c(1,1.25)), file="/home/redgar25/scratch/recluster_hepatocytes_04.jpeg", h=5,w=15)
ggsave(plot_grid(DimPlot(seu_cleaned,group.by="RNA_snn_res.0.6", label=T), DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage"), rel_widths = c(1,1.25)), file="/home/redgar25/scratch/recluster_hepatocytes_06.jpeg", h=5,w=15)

#load("/home/redgar25/scratch/tmp_hepatocytes_seurat.RData")

DotPlot(seu_cleaned, group.by = "RNA_snn_res.0.4", features = c("CYP3A4","APOE","APOA2","HAL", "RPS6"))

ggsave(DotPlot(seu_cleaned, group.by = "RNA_snn_res.0.4", features = c("CYP3A4","APOE","APOA2","HAL","RPS6","UGT2B7")), file="/home/redgar25/scratch/zonation_dot.pdf", h=5,w=10)


de<-FindMarkers(seu_cleaned, ident.1 = "6")
FeaturePlot(seu_cleaned, features = c("CYP3A4","APOE","APOA2","HAL","RPS6"))
ggsave(FeaturePlot(seu_cleaned, features = c("CYP3A4","APOE","APOA2","HAL","RPS6","percent.mt")), file="/home/redgar25/scratch/zonation_feature.jpeg", h=15,w=10)
ggsave(FeaturePlot(seu_cleaned, features = c("CYP3A4","APOE","APOA2","HAL","RPS6","MALAT1","UGT2B7","nFeature_RNA","percent.mt")), file="/home/redgar25/scratch/zonation_feature.jpeg", h=15,w=15)

de<-FindMarkers(seu_cleaned, ident.1 = "8")

save(de, file="/home/redgar25/scratch/de_cluster8.RData")

load("/home/redgar/Documents/HLiCA/data/de_cluster8.RData")

### GSEA of cluster 8 genes
source("scripts/GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("../liver_ped_map/data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")


de$gene<-rownames(de)
gene_list = de$avg_log2FC
names(gene_list) = de$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-plt_path$Enrichment

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))

ggsave(ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
         theme_bw()+ylab("")+xlab("Normalized Enrichment Score")+
         geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
         geom_hline(yintercept=16.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue")),
       file=here("figures/png/cluster8_hepatocyte_GSEA.png"), width = 15, h=10)

table(seu_cleaned$Beta.Annotation.SubLineage, seu_cleaned$RNA_snn_res.0.4)
table(seu_cleaned$STUDY, seu_cleaned$RNA_snn_res.0.4)
table(seu_cleaned$suspension_type, seu_cleaned$RNA_snn_res.0.4)
table(seu_cleaned$malat1_threshold, seu_cleaned$RNA_snn_res.0.4)
table(seu_cleaned$scmap_cell_anno, seu_cleaned$RNA_snn_res.0.4)
table(seu_cleaned$Potential.Doublets.2, seu_cleaned$RNA_snn_res.0.4)



# 
# ##############
# ## another filtering round
# ##############
# remove_second<-seu_cleaned@meta.data[which(seu_cleaned$RNA_snn_res.0.4=="6"),]
# remove_second$Gamma.Annotation<-"removed"
# save(remove_second, file="/home/redgar25/scratch/remove_second_hepatocytes_meta.RData")
# 
# 
# 
# obj <- readRDS("/home/redgar25/scratch/Hepatocyte_clean.rds")
# seu <- obj$seu
# rm(obj)
# gc()
# load("/home/redgar25/scratch/remove_second_hepatocytes_meta.RData")
# load("/home/redgar25/scratch/remove_hepatocytes_meta.RData")
# 
# cells_to_keep <- colnames(seu)[which(!(colnames(seu) %in% c(rownames(remove),rownames(remove_second))))]
# seu_cleaned<-subset(seu, cells = cells_to_keep)
# rm(seu)
# gc()
# 
# 
# dim(seu_cleaned)
# table(seu_cleaned$RNA_snn_res.0.4)
# 
# 
# ## parameters match other runs by Jordan
# seu_cleaned <- seu_cleaned %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData()
# seu_cleaned <- RunPCA(seu_cleaned, npcs = 20, verbose = FALSE)
# 
# covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')
# 
# seu_cleaned <- RunHarmony(seu_cleaned, group.by.vars = covariates, reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
# seu_cleaned <- RunUMAP(seu_cleaned, reduction = "harmony", assay = "RNA", dims = 1:20)
# seu_cleaned
# 
# seu_cleaned <- FindNeighbors(seu_cleaned, reduction = "harmony", dims = 1:20)
# seu_cleaned <- FindClusters(seu_cleaned, resolution = 0.3)
# seu_cleaned <- FindClusters(seu_cleaned, resolution = 0.4)
# seu_cleaned <- FindClusters(seu_cleaned, resolution = 0.6)
# 
# 
# ggsave(plot_grid(DimPlot(seu_cleaned,group.by="RNA_snn_res.0.3", label=T), DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage"), rel_widths = c(1,1.25)), file="/home/redgar25/scratch/recluster_hepatocytes_03_round2.jpeg", h=5,w=15)
# ggsave(plot_grid(DimPlot(seu_cleaned,group.by="RNA_snn_res.0.4", label=T), DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage"), rel_widths = c(1,1.25)), file="/home/redgar25/scratch/recluster_hepatocytes_04_round2.jpeg", h=5,w=15)
# ggsave(plot_grid(DimPlot(seu_cleaned,group.by="RNA_snn_res.0.6", label=T), DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage"), rel_widths = c(1,1.25)), file="/home/redgar25/scratch/recluster_hepatocytes_06_round2.jpeg", h=5,w=15)
# 
# 
# DotPlot(seu_cleaned, group.by = "RNA_snn_res.0.4", features = c("CYP3A4","APOE","APOA2","HAL", "RPS6"))
# 
# ggsave(DotPlot(seu_cleaned, group.by = "RNA_snn_res.0.4", features = c("CYP3A4","APOE","APOA2","HAL","RPS6")), file="/home/redgar25/scratch/zonation_dot_round2.pdf", h=5,w=10)
# ggsave(FeaturePlot(seu_cleaned, features = c("CYP3A4","APOE","APOA2","HAL","RPS6","MALAT1","UGT2B7","nFeature_RNA","percent.mt")), file="/home/redgar25/scratch/zonation_feature_round2.jpeg", h=15,w=15)


###############
## relabel
###############
unique(seu_cleaned$Beta.Annotation.SubLineage)
seu_cleaned$Gamma.Annotation<-as.character(seu_cleaned$RNA_snn_res.0.4)

seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("0","1"))]<-"Periportal Hepatocyte"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("3"))]<-"Pericentral Hepatocyte"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("2","4","5","9","11"))]<-"Ribosomal+ Hepatocyte"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("8"))]<-"UGT+ Hepatocyte"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("6"))]<-"Mito+ Hepatocyte"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("7"))]<-"SERPINE1+ Hepatocyte"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("10"))]<-"Cycling"

seu_cleaned$lineage_clusters<-seu_cleaned$RNA_snn_res.0.4

## save lineage object
save(seu_cleaned, file=here("/home/redgar25/scratch/seu_cleaned_hepatocytes.RData"))

## save meta
hepatocyte<-seu_cleaned@meta.data
save(hepatocyte, file="/home/redgar25/scratch/hepatocytes_meta.RData")



ggsave(plot_grid(DimPlot(seu_cleaned,group.by="RNA_snn_res.0.4", label=T), DimPlot(seu_cleaned, group.by = "Gamma.Annotation"), rel_widths = c(1,1.25)), file="/home/redgar25/scratch/hepatocytes_gamma.jpeg", h=5,w=15)












################## 
## jamboree feedback checks
##################
#RP emsembl IDs
RP_gene<-read.table(here("data/RP_ENSG.txt"), header=T, sep="\t")
ENSG_RP<-RP_gene$Ensembl.gene.ID[grep("^RP", RP_gene$Approved.symbol)]


## compute can

load("/home/redgar25/scratch/seu_cleaned_hepatocytes.RData")


RP<-rownames(seu_cleaned)[grep("^RP", rownames(seu_cleaned))]
noRBP <- rownames(seu_cleaned)[-(which(rownames(seu_cleaned) %in%RP))]
seu_cleaned_noRP <- subset(seu_cleaned, features = noRBP)
seu_cleaned
seu_cleaned_noRP

seu_cleaned_noRP <- seu_cleaned_noRP %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% ScaleData()
seu_cleaned_noRP <- RunPCA(seu_cleaned_noRP, npcs = 20, verbose = FALSE)

covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')

seu_cleaned_noRP <- RunHarmony(seu_cleaned_noRP, group.by.vars = covariates, reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
seu_cleaned_noRP <- RunUMAP(seu_cleaned_noRP, reduction = "harmony", assay = "RNA", dims = 1:20)
seu_cleaned_noRP

seu_cleaned_noRP <- FindNeighbors(seu_cleaned_noRP, reduction = "harmony", dims = 1:20)
seu_cleaned_noRP <- FindClusters(seu_cleaned_noRP, resolution = 0.3)
seu_cleaned_noRP <- FindClusters(seu_cleaned_noRP, resolution = 0.4)
seu_cleaned_noRP <- FindClusters(seu_cleaned_noRP, resolution = 0.6)


ggsave(plot_grid(DimPlot(seu_cleaned_noRP,group.by="RNA_snn_res.0.3", label=T), DimPlot(seu_cleaned_noRP, group.by = "Gamma.Annotation"), rel_widths = c(1,1.25)), file="/home/redgar25/scratch/recluster_hepatocytes_noRP_03.jpeg", h=5,w=15)
ggsave(plot_grid(DimPlot(seu_cleaned_noRP,group.by="RNA_snn_res.0.4", label=T), DimPlot(seu_cleaned_noRP, group.by = "Gamma.Annotation"), rel_widths = c(1,1.25)), file="/home/redgar25/scratch/recluster_hepatocytes_noRP_04.jpeg", h=5,w=15)
ggsave(plot_grid(DimPlot(seu_cleaned_noRP,group.by="RNA_snn_res.0.6", label=T), DimPlot(seu_cleaned_noRP, group.by = "Gamma.Annotation"), rel_widths = c(1,1.25)), file="/home/redgar25/scratch/recluster_hepatocytes_noRP_06.jpeg", h=5,w=15)


Idents(seu_cleaned_noRP)<-"Gamma.Annotation"
de<-FindMarkers(seu_cleaned_noRP, ident.1 = "Ribosomal+ Hepatocyte")

save(de, file="/home/redgar25/scratch/de_ribo_without_RP.RData")


## local
load("/home/redgar/Documents/HLiCA/data/de_ribo_without_RP.RData")

### GSEA of cluster 8 genes
source("scripts/GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("../liver_ped_map/data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")


de$gene<-rownames(de)
gene_list = de$avg_log2FC
names(gene_list) = de$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-plt_path$Enrichment

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=2.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))

ggsave(ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
         theme_bw()+ylab("")+xlab("Normalized Enrichment Score")+
         geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
         geom_hline(yintercept=2.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue")),
       file=here("figures/png/ribo_noRP_hepatocyte_GSEA.png"), width = 15, h=10)





load("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/hepatocytes_meta.RData")

count_plt<-as.data.frame(hepatocyte %>% dplyr::select(Gamma.Annotation,suspension_type) %>% group_by(suspension_type) %>% count(Gamma.Annotation ))

cell_count<-ggplot(count_plt, aes(Gamma.Annotation,n,  fill=suspension_type))+
  geom_bar(stat="identity", color="black")+theme_bw()+scale_fill_manual(values=c("#abdda4","#fdae61"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Cell Count")+xlab("Gamma Annotation")
cell_count
save_plts(cell_count, "suspension_effects", h=5,w=6)



#### Fancy UMAP

load(("/home/redgar25/scratch/seu_cleaned_hepatocytes.RData"))
UMAP = as.data.frame(Embeddings(seu_cleaned, reduction = "umap"))

meta<-seu_cleaned@meta.data
rm(seu_cleaned)
gc()

meta$Gamma.Annotation[which(meta$Gamma.Annotation=="Cycling Cells")]<-"Cycling"

fanciest_UMAP(meta, UMAP)
save_plts(fanciest_UMAP(meta, UMAP), "hepatocyte_fancy", w=5, h=3.5)

fanciest_UMAP(meta, UMAP, rnd_col = T)
save_plts(fanciest_UMAP(meta, UMAP,rnd_col = T), "hepatocyte_fancy_rndcol", w=5, h=3.5)


