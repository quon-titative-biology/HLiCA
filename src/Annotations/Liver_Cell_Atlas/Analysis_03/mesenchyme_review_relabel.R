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
doublets<-seu@meta.data[which(seu$seurat_clusters%in%c(5,6)),]
#doublets<-seu@meta.data[which(seu$Beta.Annotation.SubLineage=="TBC"),]
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


FeaturePlot(seu_cleaned, features = c("RELN", "NRXN1","PLP1") )
plot_grid(DimPlot(seu_cleaned, label=T),DotPlot(seu_cleaned,features = c("RELN","HHIP","COLEC10","VIPR1","LAMC3","SPON1","C7","FBLN1","RGS5","RGS6","MYH11","NRXN1","CUX2") ))
plot_grid(DimPlot(seu_cleaned, label=T),DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage",label=T),DotPlot(seu_cleaned,features = c("RELN","HHIP","COLEC10","VIPR1","LAMC3","SPON1","C7","FBLN1","RGS5","RGS6","MYH11","NRXN1","CUX2") ))

DimPlot(seu_cleaned, group.by="suspension_type",label=T)
FeaturePlot(seu_cleaned, features = c("CUX2") )



unique(seu_cleaned$Beta.Annotation.SubLineage)
seu_cleaned$Gamma.Annotation<-as.character(seu_cleaned$seurat_clusters)

seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("1","3"))]<-"Hepatic Stellate Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("0","2","5","6","8","10"))]<-"Vascular Smooth Muscle Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("4"))]<-"Portal Fibroblast"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("7"))]<-"CUX2+ Hepatic Stellate Cell"
seu_cleaned$Gamma.Annotation[which(seu_cleaned$Gamma.Annotation%in%c("9"))]<-"NRXN1+ Hepatic Stellate Cell"




DimPlot(seu_cleaned, group.by="suspension_type")
head(de)

de<-FindMarkers(seu_cleaned, ident.1 = "9")
FeaturePlot(seu_cleaned, features = c("NRXN1","NRXN3","CADM2") )

de<-FindMarkers(seu_cleaned, ident.1 = "8")
FeaturePlot(seu_cleaned, features = c("ENSG00000253295","THSD8","H2BC5") )

de<-FindMarkers(seu_cleaned, ident.1 = "10")
FeaturePlot(seu_cleaned, features = c("CYSLTR2","ADGRF5","CYGB") )


de<-FindMarkers(seu_cleaned, ident.1 = "3")
de[which(de$avg_log2FC>0),][1:10,]
FeaturePlot(seu_cleaned, features = c("HP","SERPINA1","APOC3","ALB") )
plot_grid(DimPlot(seu_cleaned, label=T),DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage",label=T),DotPlot(seu_cleaned,features = c("RELN","HHIP","COLEC10","VIPR1","LAMC3","SPON1","C7","FBLN1","RGS5","RGS6","MYH11","CUX2","HP","SERPINA1","APOC3","ALB") ))


plot_grid(DimPlot(seu_cleaned),DimPlot(seu_cleaned, group.by = "Beta.Annotation.SubLineage",label=T),DimPlot(seu_cleaned, group.by = "Gamma.Annotation",label=T), DimPlot(seu_cleaned, group.by = "RNA_snn_res.0.2",label=T))




## save mesenchyme lineage object
seu_cleaned$lineage_clusters<-seu_cleaned$RNA_snn_res.0.2
save(seu_cleaned, file=here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_mesenchyme.RData"))

## save meta
mesenchyme<-seu_cleaned@meta.data
save(mesenchyme, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/mesenchyme_meta.RData")







#####################
## schwann cells
#####################
obj <- readRDS(here("/media/redgar/Seagate Portable Drive/HLiCA/beta_annotation_objects/Mesenchyme_clean.rds"))
seu <- obj$seu
rm(obj)
gc()

DimPlot(seu, label=T)
save_plts(DimPlot(seu, label=T), "mesenchyme_clustering", h=4, w=5)

DimPlot(seu, group.by="Beta.Annotation.SubLineage")
FeaturePlot(seu, features = c("PLP1","LGI4","NRXN1") )
de<-FindMarkers(seu, ident.1 = "4")

table(seu$seurat_clusters, seu$STUDY)
table(seu$seurat_clusters, seu$suspension_type)

table(seu$STUDY)
(as.data.frame(table(seu$seurat_clusters, seu$STUDY))[which(as.data.frame(table(seu$seurat_clusters, seu$STUDY))$Var1==4),]$Freq/table(seu$STUDY))*100


### GSEA of cluster 4 genes
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
       file=here("figures/png/cluster4_mesenchyme_GSEA.png"), width = 12, h=8)




#### Fancy UMAP

load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/seu_mesenchyme.RData"))
UMAP = as.data.frame(Embeddings(seu_cleaned, reduction = "umap"))

meta<-seu_cleaned@meta.data
rm(seu_cleaned)
gc()

meta$Gamma.Annotation[which(meta$Gamma.Annotation=="Cycling Cells")]<-"Cycling"

fanciest_UMAP(meta, UMAP)
save_plts(fanciest_UMAP(meta, UMAP), "mesenchyme_fancy", w=5, h=3.5)

fanciest_UMAP(meta, UMAP, rnd_col = T)
save_plts(fanciest_UMAP(meta, UMAP,rnd_col = T), "mesenchyme_fancy_rndcol", w=5, h=3.5)

