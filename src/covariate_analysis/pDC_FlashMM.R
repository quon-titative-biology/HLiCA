### Load libraries
library(here)
library(Seurat)
library(dplyr)
library(reshape2)
library(FLASHMM)

library(ggplot2)
library(gtools)
library(cowplot)
library(RColorBrewer)
library(scales)


source("scratch/HLiCA/volcano_FlashMM.R")



covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')


##############
## Differential expression with sex 
##############

seu_cleaned<-readRDS(here("scratch/HLiCA/healthy_RNA_merged_harmonized_gamma.rds"))
seu_cleaned_pdc<-subset(seu_cleaned, subset = Gamma.Annotation == "pDC")
rm(seu_cleaned)
gc()

seu_cleaned<-seu_cleaned_pdc

## get counts
counts <- as.matrix(seu_cleaned@assays$RNA@counts)
## library size to metadata
metadata <- seu_cleaned@meta.data
metadata$libsize<-colSums(counts)

metadata$cell_id<-rownames(metadata)

## only gene expression in 10% of cells in each group
male_cells <- metadata %>% filter(donor_sex == "male") %>% pull(cell_id)
female_cells <- metadata %>% filter(donor_sex == "female") %>% pull(cell_id)

counts_male <- counts[, male_cells]
counts_female <- counts[, female_cells]

zero_percent_male <- rowSums(counts_male == 0) / ncol(counts_male)
zero_percent_female <- rowSums(counts_female == 0) / ncol(counts_female)

# Keep genes expressed in >10% of cells in both sexes
gene_to_test <- intersect(
  names(zero_percent_male)[zero_percent_male < 0.9],
  names(zero_percent_female)[zero_percent_female < 0.9]
)

counts<-counts[gene_to_test,]

metadata$donor_sex<-as.factor(metadata$donor_sex)
metadata$donor_sex <- relevel(metadata$donor_sex, ref = "female")



c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')
# 'donor_uuid', 'library_alias' as random effects?

##Gene expression matrix, Y = log2(1 + counts)
Y <- log2(1 + counts)

##Design matrix for fixed effects
#X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = meta)
#X <- model.matrix(~ 0 + log(libsize) + donor_sex + donor_sex:Gamma.Annotation +  STUDY + suspension_type + assay + donor_uuid , data = metadata)
X <- model.matrix(~ log(libsize) + donor_sex + STUDY + suspension_type, data = metadata)

##Design matrix for random effects
#Z <- model.matrix(~ 0 + as.factor(sam), data = dat$meta)
Z <- model.matrix(~ 0 + as.factor(library_alias), data = metadata)

##Dimension of random effects
d <- ncol(Z)


max.iter <- 100
epsilon <- 1e-5
fit <- lmmfit(Y, X, Z, d = d, max.iter = max.iter, epsilon = epsilon)

gc <- (apply(abs(fit$dlogL), 2, max) < epsilon) 
sum(gc) 



## contrast for sex (summed across all cell types)
contrast <- cbind("FvsM" = numeric(nrow(fit$coef)))
index <- grep("donor_sexmale", rownames(fit$coef))
contrast[index, ] <- 1/length(index)

##Test the contrast.
test <- lmmtest(fit, contrast = contrast)
head(test)

test<-as.data.frame(test)

out <- data.frame(
  gene = rownames(test), 
  coef = test$FvsM_coef, p = test$FvsM_p)


##Adjust p-values by FDR.
out$FDR <- p.adjust(out$p, method = "fdr")
# flip direction for volcano
out$coef <- (-out$coef)


##The DE genes with FDR < 0.05 and abs(logFC) > 0.5
out <- out[order(out$p), ]
rownames(out) <- NULL
out<-out[which(!(is.na(out$p))),]
sig_de<-out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]
nrow(sig_de)


### pathway for all
source("scratch/HLiCA/00_GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("scratch/HLiCA/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")


gene_list = out$coef
names(gene_list) = out$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
res<-GSEA(gene_list, GO_file, pval = 0.05)


plt_path<-res$Results

plt_path<-rbind(plt_path[1:15,], plt_path[order(plt_path$NES),][1:15,])
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
plt_path$Enrichment<-as.factor(plt_path$Enrichment)
levels(plt_path$Enrichment)<-c("Enriched in Males", "Enriched in Females")

plt_path_save<-plt_path[,c("pathway","Enrichment","pval","padj","ES","NES","nMoreExtreme","size","leadingEdge")]
plt_path_save$leadingEdge<-sapply(1:nrow(plt_path_save), function(x) paste0(plt_path_save$leadingEdge[x][[1]], collapse = ", "))
write.csv(plt_path_save, file=paste("scratch/HLiCA/age_DEG_GSEA/","pDC", "_sex_GSEA.csv", sep=""), row.names = FALSE)

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

pdc_pathway<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment), shape=21)+
  theme_bw()+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=15.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))+
  ggtitle("pDC")
pdc_pathway
ggsave(pdc_pathway, file=paste("scratch/HLiCA/","pDC","_sex_GSEA.png",sep=""), w=14,h=7)


pdc_volcano<-makeVolcano(out, 0.5, 0.05, paste("Differential expression between sexes\n","pDC", sep=""), round(max(abs(range(out$coef))),1)+0.2)
pdc_volcano
ggsave(pdc_volcano, file=paste("scratch/HLiCA/","pDC","_sex_volcano.png",sep=""), w=10,h=7)

sig_de<-out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]
write.csv(sig_de, file=paste("scratch/HLiCA/age_DEG_GSEA/","pDC", "_sex_DEG.csv", sep=""), row.names = FALSE)





##########################################################################################################################
## Differential expression with AGE all pdc
##########################################################################################################################


## age catagories
age<-read.csv("scratch/HLiCA/HLiCA_scadmix_geneticAncestry_plus_age_catagories_manual_with_library_uuid.csv")

meta_add<-seu_cleaned@meta.data[,c("library_uuid","orig.ident")]
meta_add$index<-rownames(meta_add)
meta_add<-merge(meta_add, age, by="orig.ident")
dim(meta_add)

meta_add<-meta_add[match(colnames(seu_cleaned), meta_add$index),]
identical(colnames(seu_cleaned), meta_add$index)

seu_cleaned <- AddMetaData(seu_cleaned, metadata = meta_add)

#### Age continuum

seu_cleaned$age_con<-as.factor(seu_cleaned$donor_age_cat)
levels(seu_cleaned$age_con)<-c( "<19", NA, "19-30","31-40","41-50","51-60","61-70","71-80","81-90" )
seu_cleaned$age_con<-factor(seu_cleaned$age_con, levels = c(  NA,"<19",  "19-30"   ,  "31-40"   ,  "41-50"   ,  "51-60" ,    "61-70"   ,  "71-80"   ,  "81-90" ))
seu_cleaned$age_int<-as.numeric(seu_cleaned$age_con)

table(seu_cleaned$age_int, seu_cleaned$age_con)


with_age<-colnames(seu_cleaned)[which(!(is.na(seu_cleaned$age_int)))]
seu_cleaned$cell<-colnames(seu_cleaned)
seu_cleaned_age<-subset(seu_cleaned, subset = cell %in% with_age)


counts <- as.matrix(seu_cleaned_age@assays$RNA@counts)

metadata <- seu_cleaned_age@meta.data
metadata$libsize<-colSums(counts)

metadata$age_int_scaled <- scale(metadata$age_int, center = TRUE, scale = FALSE)

## only gene expression in 10% of cells in each group
zero_percent <- rowSums(counts == 0) / ncol(counts)

# Keep genes expressed in >10% of cells in both sexes
gene_to_test <- names(zero_percent)[zero_percent < 0.95]

counts<-counts[gene_to_test,]


metadata %>%
  select(age_int_scaled, library_alias) %>%
  distinct() %>%
  count(age_int_scaled, name = "n_unique_libraries")
# 
# 
# c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')
# # 'donor_uuid', 'library_alias' as random effects?
# 
# ##Gene expression matrix, Y = log2(1 + counts)
# Y <- log2(1 + counts)
# 
# ##Design matrix for fixed effects
# #X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = meta)
# #X <- model.matrix(~ 0 + log(libsize) + age_int, data = metadata)
# #X <- model.matrix(~ 0 + log(libsize) + age_int_scaled, data = metadata)
# X <- model.matrix(~ 0 + log(libsize) + age_int_scaled, data = metadata)
# 
# 
# ##Design matrix for random effects
# #Z <- model.matrix(~ 0 + as.factor(sam), data = dat$meta)
# Z <- model.matrix(~ 0 + as.factor(library_alias), data = metadata)
# 
# ##Dimension of random effects
# d <- ncol(Z)
# 
# 
# max.iter <- 100
# epsilon <- 1e-5
# fit <- lmmfit(Y, X, Z, d = d, max.iter = max.iter, epsilon = epsilon)
# 
# gc <- (apply(abs(fit$dlogL), 2, max) < epsilon) 
# sum(gc) 
# 
# 
# 
# ##############
# ## Any cell type specific differences?
# ##############
# celltypes<- unique(seu_cleaned$Gamma.Annotation)
# 
# age_out_celltype<-lapply(celltypes, function(celltype){
#   out <- data.frame(
#     gene = colnames(fit$p), 
#     coef = fit$coef[which(rownames(fit$coef) == paste("Gamma.Annotation",celltype,":age_int_scaled",sep="")),], p = fit$p[which(rownames(fit$p) == paste("Gamma.Annotation",celltype,":age_int_scaled",sep="")),])
#   
#   ##Adjust p-values by FDR.
#   out$FDR <- p.adjust(out$p, method = "fdr")
#   
#   ##The DE genes with FDR < 0.05 and abs(logFC) > 0.5
#   out <- out[order(out$p), ]
#   rownames(out) <- NULL
#   out
# })
# 
# names(age_out_celltype)<-celltypes
# 
# 
# 
# ########### significant gene count
# do.call(rbind,lapply(1:length(age_out_celltype), function(x){
#   age_out_celltype[[x]]
#   sig_de<-age_out_celltype[[x]][(age_out_celltype[[x]]$FDR < 0.05) & (abs(age_out_celltype[[x]]$coef) > 0.5) , ]
#   data.frame(celltype=names(age_out_celltype)[[x]], number_sig=nrow(sig_de))
# }))
# 
# 
# ### pathway for all
# source("scratch/HLiCA/00_GSEA_function.R")
# #http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
# GO_file = here("scratch/HLiCA/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")
# 
# 
# GSEA_age_celltype<-lapply(1:length(age_out_celltype), function(x){
#   print(names(age_out_celltype)[x])
#   out<-age_out_celltype[[x]]
#   gene_list = out$coef
#   names(gene_list) = out$gene
#   gene_list = sort(gene_list, decreasing = TRUE)
#   gene_list = gene_list[!duplicated(names(gene_list))]
#   GSEA(gene_list, GO_file, pval = 0.05)
# })
# 
# 
# 
# 
# 
# #### plot paths if anything significant
# save_DEG_plots_age<-function(celltype_num){
#   print(names(age_out_celltype)[celltype_num])
#   
#   plt_path<-GSEA_age_celltype[[celltype_num]]$Results
#   if(nrow(plt_path)==0){}else{
#     plt_path<-rbind(plt_path[1:15,], plt_path[order(plt_path$NES),][1:15,])
#     plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
#     plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
#     plt_path$Enrichment<-as.factor(plt_path$Enrichment)
#     levels(plt_path$Enrichment)<-c("Enriched in genes decreasing with age", "Enriched in genes increasing with age")
#     
#     plt_path_save<-plt_path[,c("pathway","Enrichment","pval","padj","ES","NES","nMoreExtreme","size","leadingEdge")]
#     plt_path_save$leadingEdge<-sapply(1:nrow(plt_path_save), function(x) paste0(plt_path_save$leadingEdge[x][[1]], collapse = ", "))
#     write.csv(plt_path_save, file=paste("scratch/HLiCA/age_DEG_GSEA/",names(age_out_celltype)[celltype_num], "_age_GSEA.csv", sep=""), row.names = FALSE)
#     
#     plt_path$direction_label<-as.factor(plt_path$Enrichment)
#     levels(plt_path$direction_label)<-c(0.1,-0.1)
#     plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
#     
#     pdc_pathway<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment), shape=21)+
#       theme_bw()+ylab("")+xlab("Normalized Enrichment Score")+
#       geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
#       geom_hline(yintercept=15.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))+
#       ggtitle(names(age_out_celltype)[celltype_num])
#     pdc_pathway
#     ggsave(pdc_pathway, file=paste("scratch/HLiCA/",names(age_out_celltype)[celltype_num],"_age_GSEA.png",sep=""), w=14,h=7)
#   }
#   
#   pdc_volcano<-makeVolcano_age(age_out_celltype[[celltype_num]], 0.5, 0.05, paste("Differential expression with age\n",names(age_out_celltype)[celltype_num], sep=""), round(max(abs(range(age_out_celltype[[celltype_num]]$coef))),1)+0.2)
#   pdc_volcano
#   ggsave(pdc_volcano, file=paste("scratch/HLiCA/",names(age_out_celltype)[celltype_num],"_age_volcano.png",sep=""), w=10,h=7)
#   
#   sig_de<-age_out_celltype[[celltype_num]][(age_out_celltype[[celltype_num]]$FDR < 0.05) & (abs(age_out_celltype[[celltype_num]]$coef) > 0.5) , ]
#   write.csv(sig_de, file=paste("scratch/HLiCA/age_DEG_GSEA/",names(age_out_celltype)[celltype_num], "_age_DEG.csv", sep=""), row.names = FALSE)
# }
# 
# lapply(1:length(celltypes), function(x) save_DEG_plots_age(x))
# 
# 
# 
# 
# 
# # ### smooth line violin plot
# # 
# # data_plot<-FetchData(seu_cleaned_age, vars = sig_de$gene)
# # data_plot$cell_id<-rownames(data_plot)
# # data_plot<-melt(data_plot)
# # meta<-seu_cleaned_age@meta.data
# # meta$cell_id<-rownames(meta)
# # 
# # data_plot<-merge(data_plot, meta, by="cell_id")
# # 
# # 
# # Means <- data_plot %>% group_by(age_con, variable) %>% 
# #   summarize(Avg = mean(value))
# # 
# # age_lymphocyte<-ggplot() + 
# #   geom_violin(data = data_plot, mapping = aes(x = age_con, y = value), color="grey", fill="grey") +
# #   geom_point(data = Means, mapping = aes(x = age_con, y = Avg), color="red") +
# #   geom_line(data = Means, mapping = aes(x = age_con, y = Avg, group=variable), color="red")+
# #   geom_smooth(data = Means, mapping = aes(x = age_con, y = Avg, group=variable), method="lm")+
# #   facet_wrap(~variable, scales = "free_y")+theme_bw()
# # age_lymphocyte
# # 
# # ggsave(age_lymphocyte, file="scratch/HLiCA/age_lymphocyte_genes.png", h=10, w=14)
# # 
# # 
# # 
# # data_plot<-FetchData(seu_cleaned_age, vars = sig_de$gene)
# # data_plot$cell_id<-rownames(data_plot)
# # data_plot<-melt(data_plot)
# # meta<-seu_cleaned_age@meta.data
# # meta$cell_id<-rownames(meta)
# # 
# # data_plot<-merge(data_plot, meta, by="cell_id")
# # 
# # 
# # Means <- data_plot %>% group_by(age_con, variable) %>% 
# #   summarize(Avg = mean(value))
# # 
# # age_lymphocyte<-ggplot() + 
# #   geom_violin(data = data_plot, mapping = aes(x = age_con, y = value), color="grey", fill="grey") +
# #   geom_point(data = Means, mapping = aes(x = age_con, y = Avg), color="red") +
# #   geom_line(data = Means, mapping = aes(x = age_con, y = Avg, group=variable), color="red")+
# #   geom_smooth(data = Means, mapping = aes(x = age_con, y = Avg, group=variable), method="lm")+
# #   facet_wrap(~variable, scales = "free_y")+theme_bw()
# # age_lymphocyte
# # 
# # ggsave(age_lymphocyte, file="scratch/HLiCA/age_lymphocyte_genes.png", h=10, w=14)
# # 
# # 
