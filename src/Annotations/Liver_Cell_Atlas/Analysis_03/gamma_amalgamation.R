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
library(data.table)


source("scripts/00_pretty_plots.R")


load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/myeloid_meta.RData"))
load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/cholangiocyte_meta.RData"))
load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/endo_meta.RData"))
load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/mesenchyme_meta.RData"))
load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/hepatocytes_meta.RData"))
load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/lymphocyte_meta.RData"))
load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/pDC_meta.RData"))

load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/doublets_cholangiocytes_meta.RData"))
doublets_cho<-doublets
load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/doublets_lymphocytes_meta.RData"))
doublets_lymp<-doublets
load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/doublets_mesenchyme_meta.RData"))
doublets_mes<-doublets
load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/doublets_myeloid_meta.RData"))
doublets_mye<-doublets

load(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/remove_hepatocytes_meta.RData"))
removed_hepatocytes<-remove


col_keep<-c(  "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt" ,      "cell_barcode" ,      "S.Score" ,   "G2M.Score" ,"Phase",
              "SAMPLE","STUDY", "tissue" ,"suspension_type"  ,   "assay", "sequencing_run","library_alias",
              "library_uuid" ,"donor_uuid" ,"donor_sex", "donor_age" ,"donor_ethnicity" ,
              "condition", "malat1_threshold" ,"Alpha.Annotations" ,"Beta.Annotation.SubLineage" ,"Gamma.Annotation")

meta_HLiCA<-rbind(myeloid[,col_keep],
                  cholangiocyte[,col_keep],
                  endo[,col_keep],
                  mesenchyme[,col_keep],
                  lymphocyte[,col_keep],
                  hepatocyte[,col_keep],
                  pDC[,col_keep])


unique(meta_HLiCA$Gamma.Annotation)
save(meta_HLiCA, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/gamma_amalgamation.RData")


meta_HLiCA_withRemoved<-rbind(myeloid[,col_keep],
                  cholangiocyte[,col_keep],
                  endo[,col_keep],
                  mesenchyme[,col_keep],
                  lymphocyte[,col_keep],
                  hepatocyte[,col_keep],
                  pDC[,col_keep],
                  doublets_cho[,col_keep],
                  doublets_lymp[,col_keep],
                  doublets_mes[,col_keep],
                  doublets_mye[,col_keep],
                  removed_hepatocytes[,col_keep])

save(meta_HLiCA_withRemoved, file="/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/gamma_amalgamation_withRemoved.RData")






#################
## load original object
#################
seu <- readRDS("/home/redgar25/scratch/healthy_RNA_merged_harmonized.rds") # Load data
load("/home/redgar25/scratch/gamma_amalgamation.RData")

meta_HLiCA$Gamma.Annotation[which(meta_HLiCA$Gamma.Annotation=="Cycling Cells")]<-"Cycling"


#################
## remove cells filtered in gamma
#################
seu

cells_to_keep <- colnames(seu)[which(colnames(seu) %in% rownames(meta_HLiCA))]
length(cells_to_keep)
seu_cleaned<-subset(seu, cells = cells_to_keep)
seu_cleaned
rm(seu)
gc()

meta_HLiCA<-meta_HLiCA[match(colnames(seu_cleaned), rownames(meta_HLiCA)),]
identical(colnames(seu_cleaned), rownames(meta_HLiCA))

seu_cleaned <- AddMetaData(seu_cleaned, metadata = meta_HLiCA)
seu_cleaned@meta.data<-seu_cleaned@meta.data[,col_keep]

covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')

seu_cleaned <- seu_cleaned %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% ScaleData()
seu_cleaned <- RunPCA(seu_cleaned, assay="RNA", features=rownames(seu_cleaned), npcs=40)
covariates
seu_harmonized <- RunHarmony(seu_cleaned, group.by.vars = covariates, reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
seu_harmonized <- RunUMAP(seu_harmonized, reduction = "harmony", assay = "RNA", dims = 1:40)
seu_harmonized

HARMONY = Embeddings(seu_harmonized, reduction = "harmony")
UMAP = Embeddings(seu_harmonized, reduction = "umap")
META = seu_harmonized@meta.data
fwrite(as.data.frame(HARMONY), file.path('/home/redgar25/scratch/harmony_embedding_gamma.csv'), sep = ',', quote=F, row.names = T, col.names = T)
fwrite(as.data.frame(UMAP), file.path('/home/redgar25/scratch/harmony_umap_gamma.csv'), sep = ',', quote=F, row.names = T, col.names = T)
fwrite(as.data.frame(META), file.path('/home/redgar25/scratch/healthy_metadata_gamma.csv'), sep = ',', quote=F, row.names = T, col.names = T)
saveRDS(seu_harmonized, file.path('/home/redgar25/scratch/healthy_RNA_merged_harmonized_gamma.rds'))



#################
## with doublets for extra analysis (also the ones Jordan removed?) - how to mark them systematically
#################
seu <- readRDS("/home/redgar25/scratch/healthy_RNA_merged_harmonized.rds") # Load data
load(here("/home/redgar25/scratch/annotation_tracking.RData"))
load("/home/redgar25/scratch/gamma_amalgamation.RData")



kept<-meta_HLiCA[,c("Alpha.Annotations","Beta.Annotation.SubLineage","Gamma.Annotation")]
kept$filtered<-NA

removed<-df_forsankey[which(!(df_forsankey$X%in%rownames(meta_HLiCA))),]
rownames(removed)<-removed$X
removed$filtered<-"Clustered by technical factors once sorted into alpha lineage"
removed$filtered[which(removed$Alpha.Annotations.Doublets=="Potential Doublet")]<-"Doublet seen in alpha annotation"
removed$filtered[which(removed$Beta.Annotations=="TBC")]<-"Doublet seen in beta annotation"
removed$filtered[which(removed$Alpha.Annotations=="Erythrocyte")]<-"Erythrocyte"
removed$filtered[which(removed$Alpha.Annotations=="Doublet")]<-"Doublet seen in alpha annotation"
removed$Gamma.Annotation<-NA
removed<-removed[,c("Alpha.Annotations","Beta.Annotation.SubLineage","Gamma.Annotation","filtered")]


meta_raw<-rbind(kept, removed)

meta_raw<-meta_raw[match(colnames(seu), rownames(meta_raw)),]
identical(colnames(seu), rownames(meta_raw))

meta_raw$Gamma.Annotation[which(meta_raw$Gamma.Annotation=="Cycling Cells")]<-"Cycling"

seu <- AddMetaData(seu, metadata = meta_raw)
seu@meta.data<-seu@meta.data[,which(colnames(seu@meta.data)%in%c(col_keep,"filtered"))]

seu


HARMONY = Embeddings(seu, reduction = "harmony")
UMAP = Embeddings(seu, reduction = "umap")
META = seu@meta.data
fwrite(as.data.frame(HARMONY), file.path('/home/redgar25/scratch/harmony_embedding_gamma_withRemoved.csv'), sep = ',', quote=F, row.names = T, col.names = T)
fwrite(as.data.frame(UMAP), file.path('/home/redgar25/scratch/harmony_umap_gamma_withRemoved.csv'), sep = ',', quote=F, row.names = T, col.names = T)
fwrite(as.data.frame(META), file.path('/home/redgar25/scratch/healthy_metadata_gamma_withRemoved.csv'), sep = ',', quote=F, row.names = T, col.names = T)
saveRDS(seu, file.path('/home/redgar25/scratch/healthy_RNA_merged_harmonized_gamma_withRemoved.rds'))

