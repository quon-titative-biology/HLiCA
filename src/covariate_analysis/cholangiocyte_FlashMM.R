### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(reshape2)
library(gtools)
library(cowplot)
library(RColorBrewer)
library(FLASHMM)
library(ggsignif)


source("scratch/HLiCA/volcano_FlashMM.R")



covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')


##############
## Differential expression with sex 
##############

load(here("scratch/HLiCA/seu_cholangiocyte.RData"))


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

metadata$donor_sex<-as.character(metadata$donor_sex)

c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')
# 'donor_uuid', 'library_alias' as random effects?

##Gene expression matrix, Y = log2(1 + counts)
Y <- log2(1 + counts)

##Design matrix for fixed effects
#X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = meta)
#X <- model.matrix(~ 0 + log(libsize) + donor_sex + donor_sex:Gamma.Annotation +  STUDY + suspension_type + assay + donor_uuid , data = metadata)
X <- model.matrix(~ 0 + log(libsize) + Gamma.Annotation + Gamma.Annotation:donor_sex +  STUDY + suspension_type, data = metadata)

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
index <- grep(":", rownames(fit$coef))
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
out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]

sig_de<-out[(out$FDR < 0.05) & (abs(out$coef) > 0.25) , ]
nrow(sig_de)


cholangiocyte_volcano<-makeVolcano(out, 0.25, 0.05, "Differential expression between sexes", 0.6)
cholangiocyte_volcano
ggsave(cholangiocyte_volcano, file="scratch/HLiCA/cholangiocyte_sex_volcano.png", w=10,h=7)


VlnPlot(seu_cleaned, features = c("KDM6A","TTTY14"), log = T, group.by = "donor_sex")
VlnPlot(seu_cleaned, features = sig_de$gene[1:10], group.by = "donor_sex")


source("scratch/HLiCA/00_GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("scratch/HLiCA/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

### sex
gene_list = out$coef
names(gene_list) = out$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
res = GSEA(gene_list, GO_file, pval = 0.05)

#### plot paths is anything significant

plt_path<-res$Results

### no sig pathways


##############
## Any cell type specific differences?
##############
celltypes<- unique(seu_cleaned$Gamma.Annotation)

sex_out_celltype<-lapply(celltypes, function(celltype){
  out <- data.frame(
    gene = colnames(fit$p), 
    coef = fit$coef[which(rownames(fit$coef) == paste("Gamma.Annotation",celltype,":donor_sexmale",sep="")),], p = fit$p[which(rownames(fit$p) == paste("Gamma.Annotation",celltype,":donor_sexmale",sep="")),])
  
  ##Adjust p-values by FDR.
  out$FDR <- p.adjust(out$p, method = "fdr")
  # flip direction for volcano
  out$coef <- (-out$coef)
  
  ##The DE genes with FDR < 0.05 and abs(logFC) > 0.5
  out <- out[order(out$p), ]
  rownames(out) <- NULL
  out
})

names(sex_out_celltype)<-celltypes


head(sex_out_celltype[["Large Mucus Secreting"]])

sig_de<-out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]
nrow(sig_de)


cholangiocyte_volcano<-makeVolcano(sex_out_celltype[["Large Mucus Secreting"]], 0.5, 0.05, "Differential expression between sexes", round(abs(max(range(sex_out_celltype[["Large Mucus Secreting"]]$coef))),1))
cholangiocyte_volcano
#ggsave(cholangiocyte_volcano, file="scratch/HLiCA/cholangiocyte_sex_volcano.png", w=10,h=7)


VlnPlot(seu_cleaned, features = c("TTTY14","FLRT2"), log = T, split.by = "donor_sex",  group.by = "Gamma.Annotation")
VlnPlot(seu_cleaned, features = sig_de$gene[1:10], group.by = "donor_sex")

########### significant gene count
do.call(rbind,lapply(1:length(sex_out_celltype), function(x){
  sex_out_celltype[[x]]
  sig_de<-sex_out_celltype[[x]][(sex_out_celltype[[x]]$FDR < 0.05) & (abs(sex_out_celltype[[x]]$coef) > 0.5) , ]
  data.frame(celltype=names(sex_out_celltype)[[x]], number_sig=nrow(sig_de))
}))




### pathway for all
source("scratch/HLiCA/00_GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("scratch/HLiCA/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")


GSEA_sex_celltype<-lapply(1:length(sex_out_celltype), function(x){
  print(names(sex_out_celltype)[x])
  out<-sex_out_celltype[[x]]
  gene_list = out$coef
  names(gene_list) = out$gene
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  GSEA(gene_list, GO_file, pval = 0.05)
})





#### plot paths if anything significant
# sig paths 2, 4, 5

save_DEG_plots<-function(celltype_num){
  print(names(sex_out_celltype)[celltype_num])
  
  plt_path<-GSEA_sex_celltype[[celltype_num]]$Results
  if(nrow(plt_path)==0){}else{
    plt_path<-rbind(plt_path[1:15,], plt_path[order(plt_path$NES),][1:15,])
    plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
    plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
    plt_path$Enrichment<-as.factor(plt_path$Enrichment)
    levels(plt_path$Enrichment)<-c("Enriched in Males", "Enriched in Females")
    
    plt_path_save<-plt_path[,c("pathway","Enrichment","pval","padj","ES","NES","nMoreExtreme","size","leadingEdge")]
    plt_path_save$leadingEdge<-sapply(1:nrow(plt_path_save), function(x) paste0(plt_path_save$leadingEdge[x][[1]], collapse = ", "))
    write.csv(plt_path_save, file=paste("scratch/HLiCA/age_DEG_GSEA/",names(age_out_celltype)[celltype_num], "_sex_GSEA.csv", sep=""), row.names = FALSE)
    
    plt_path$direction_label<-as.factor(plt_path$Enrichment)
    levels(plt_path$direction_label)<-c(0.1,-0.1)
    plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
    
    cholangiocyte_pathway<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment), shape=21)+
      theme_bw()+ylab("")+xlab("Normalized Enrichment Score")+
      geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
      geom_hline(yintercept=15.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))+
      ggtitle(names(sex_out_celltype)[celltype_num])
    cholangiocyte_pathway
    ggsave(cholangiocyte_pathway, file=paste("scratch/HLiCA/",names(sex_out_celltype)[celltype_num],"_sex_GSEA.png",sep=""), w=14,h=7)
  }
  
  cholangiocyte_volcano<-makeVolcano(sex_out_celltype[[celltype_num]], 0.5, 0.05, paste("Differential expression between sexes\n",names(sex_out_celltype)[celltype_num], sep=""), round(max(abs(range(sex_out_celltype[[celltype_num]]$coef))),1)+0.2)
  cholangiocyte_volcano
  ggsave(cholangiocyte_volcano, file=paste("scratch/HLiCA/",names(sex_out_celltype)[celltype_num],"_sex_volcano.png",sep=""), w=10,h=7)
  
  sig_de<-sex_out_celltype[[celltype_num]][(sex_out_celltype[[celltype_num]]$FDR < 0.05) & (abs(sex_out_celltype[[celltype_num]]$coef) > 0.5) , ]
  write.csv(sig_de, file=paste("scratch/HLiCA/age_DEG_GSEA/",names(sex_out_celltype)[celltype_num], "_sex_DEG.csv", sep=""), row.names = FALSE)
}


lapply(1:length(celltypes), function(x) save_DEG_plots(x))











##########################################################################################################################
## Differential expression with AGE all cholangiocytes
##########################################################################################################################

load(here("scratch/HLiCA/seu_cholangiocyte.RData"))


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


c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')
# 'donor_uuid', 'library_alias' as random effects?

##Gene expression matrix, Y = log2(1 + counts)
Y <- log2(1 + counts)

##Design matrix for fixed effects
#X <- model.matrix(~ 0 + log(libsize) + cls + cls:trt, data = meta)
#X <- model.matrix(~ 0 + log(libsize) + age_int, data = metadata)
#X <- model.matrix(~ 0 + log(libsize) + age_int_scaled, data = metadata)
X <- model.matrix(~ 0 + log(libsize) + Gamma.Annotation + Gamma.Annotation:age_int_scaled, data = metadata)


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


# 
# ##(3) The DE genes specific to a cell-type
# ##Coefficients, t-values, and p-values for the genes specific to a cell-type.
# index <- grep("age_int", rownames(fit$coef))
# ce <- fit$coef[index, gc]
# tv <- fit$t[index, gc]
# pv <- fit$p[index, gc]
# 
# out <- data.frame(
#   gene = names(ce), 
#   coef = ce , t = tv, p = pv)
# 
# hist(out$coef)
# 
# out$delta_age<-sapply(1:nrow(out), function(x) (out$coef[x]*(max(metadata$age_int_scaled))) - (out$coef[x]*(min(metadata$age_int_scaled))))
# 
# 
# 
# 
# 
# ##Adjust p-values by FDR.
# out$FDR <- p.adjust(out$p, method = "fdr")
# 
# ##The DE genes with FDR < 0.05 and abs(logFC) > 0.5
# out <- out[order(out$p), ]
# rownames(out) <- NULL
# out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]
# 
# sig_de<-out[(out$FDR < 0.05) & (abs(out$coef) > 0.5) , ]
# nrow(sig_de)
# 
# stat_sig<-out[(out$FDR < 0.05),]
# 
# sig_de<-out[(out$FDR < 0.05) & (abs(out$delta_age) > 1) , ]
# nrow(sig_de)
# 
# 
# head(sig_de[order(sig_de$coef),])
# 
# VlnPlot(seu_cleaned_age, features = "SAA1", group.by = "age_int" )
# VlnPlot(seu_cleaned_age, features = "KLF12", group.by = "age_int" )


##############
## Any cell type specific differences?
##############
celltypes<- unique(seu_cleaned$Gamma.Annotation)

age_out_celltype<-lapply(celltypes, function(celltype){
  out <- data.frame(
    gene = colnames(fit$p), 
    coef = fit$coef[which(rownames(fit$coef) == paste("Gamma.Annotation",celltype,":age_int_scaled",sep="")),], p = fit$p[which(rownames(fit$p) == paste("Gamma.Annotation",celltype,":age_int_scaled",sep="")),])
  
  ##Adjust p-values by FDR.
  out$FDR <- p.adjust(out$p, method = "fdr")
  
  ##The DE genes with FDR < 0.05 and abs(logFC) > 0.5
  out <- out[order(out$p), ]
  rownames(out) <- NULL
  out
})

names(age_out_celltype)<-celltypes

head(age_out_celltype[["Large Mucus Secreting"]])



cholangiocyte_volcano<-makeVolcano_age(age_out_celltype[["Large Mucus Secreting"]], 0.5, 0.05, "Differential expression between sexes", round(abs(max(range(age_out_celltype[["Large Mucus Secreting"]]$coef))),1))
cholangiocyte_volcano
#ggsave(cholangiocyte_volcano, file="scratch/HLiCA/cholangiocyte_sex_volcano.png", w=10,h=7)

VlnPlot(seu_cleaned, features = c("AGR2"), log = T, split.by = "age_con",  group.by = "Gamma.Annotation")+scale_fill_manual(values=c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84', "grey",   "grey40" ,  "#74878a"))
ggsave("scratch/HLiCA/AGR2_mucus_screting.png", w=10, h=5)


########### significant gene count
do.call(rbind,lapply(1:length(age_out_celltype), function(x){
  age_out_celltype[[x]]
  sig_de<-age_out_celltype[[x]][(age_out_celltype[[x]]$FDR < 0.05) & (abs(age_out_celltype[[x]]$coef) > 0.5) , ]
  data.frame(celltype=names(age_out_celltype)[[x]], number_sig=nrow(sig_de))
}))



### pathway for all
source("scratch/HLiCA/00_GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("scratch/HLiCA/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")


GSEA_age_celltype<-lapply(1:length(age_out_celltype), function(x){
  print(names(age_out_celltype)[x])
  out<-age_out_celltype[[x]]
  gene_list = out$coef
  names(gene_list) = out$gene
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  GSEA(gene_list, GO_file, pval = 0.05)
})





#### plot paths if anything significant
# sig paths 2, 4, 5

save_DEG_plots_age<-function(celltype_num){
  print(names(age_out_celltype)[celltype_num])
  
  plt_path<-GSEA_age_celltype[[celltype_num]]$Results
  plt_path<-rbind(plt_path[1:15,], plt_path[order(plt_path$NES),][1:15,])
  plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
  plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
  plt_path$Enrichment<-as.factor(plt_path$Enrichment)
  levels(plt_path$Enrichment)<-c("Enriched in genes decreasing with age", "Enriched in genes increasing with age")
  
  plt_path_save<-plt_path[,c("pathway","Enrichment","pval","padj","ES","NES","nMoreExtreme","size","leadingEdge")]
  plt_path_save$leadingEdge<-sapply(1:nrow(plt_path_save), function(x) paste0(plt_path_save$leadingEdge[x][[1]], collapse = ", "))
  write.csv(plt_path_save, file=paste("scratch/HLiCA/age_DEG_GSEA/",names(age_out_celltype)[celltype_num], "_age_GSEA.csv", sep=""), row.names = FALSE)
  
  plt_path$direction_label<-as.factor(plt_path$Enrichment)
  levels(plt_path$direction_label)<-c(0.1,-0.1)
  plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
  
  cholangiocyte_pathway<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment), shape=21)+
    theme_bw()+ylab("")+xlab("Normalized Enrichment Score")+
    geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
    geom_hline(yintercept=15.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))+
    ggtitle(names(age_out_celltype)[celltype_num])
  cholangiocyte_pathway
  ggsave(cholangiocyte_pathway, file=paste("scratch/HLiCA/",names(age_out_celltype)[celltype_num],"_age_GSEA.png",sep=""), w=14,h=7)
  
  cholangiocyte_volcano<-makeVolcano_age(age_out_celltype[[celltype_num]], 0.5, 0.05, paste("Differential expression between sexes\n",names(age_out_celltype)[celltype_num], sep=""), round(max(abs(range(age_out_celltype[[celltype_num]]$coef))),1)+0.2)
  cholangiocyte_volcano
  ggsave(cholangiocyte_volcano, file=paste("scratch/HLiCA/",names(age_out_celltype)[celltype_num],"_age_volcano.png",sep=""), w=10,h=7)
  
  sig_de<-age_out_celltype[[celltype_num]][(age_out_celltype[[celltype_num]]$FDR < 0.05) & (abs(age_out_celltype[[celltype_num]]$coef) > 0.5) , ]
  write.csv(sig_de, file=paste("scratch/HLiCA/age_DEG_GSEA/",names(age_out_celltype)[celltype_num], "_age_DEG.csv", sep=""), row.names = FALSE)
}

save_DEG_plots_age(1)
save_DEG_plots_age(2)
save_DEG_plots_age(3)
save_DEG_plots_age(4)
save_DEG_plots_age(5)

###############
## Sig genes all cell types
###############
sig_genes_relaxed<-do.call(rbind,lapply(1:length(age_out_celltype), function(celltype_num){
  print(names(age_out_celltype)[celltype_num])
  sig<-age_out_celltype[[celltype_num]][(age_out_celltype[[celltype_num]]$FDR < 0.05) & 
                                          (abs(age_out_celltype[[celltype_num]]$coef) > 0.1) , ]
  if(nrow(sig)>0){
    sig$Gamma.Annotation<-names(age_out_celltype)[celltype_num]
    sig}
}))


mucus_sig<-age_out_celltype[[5]][(age_out_celltype[[5]]$FDR < 0.05) & (abs(age_out_celltype[[5]]$coef) > 0.25) , ]
mucus_sig[rev(order(mucus_sig$coef)),]
age_out_celltype[[5]][grep("EZR",age_out_celltype[[5]]$gene),]

age_out_celltype[[5]][which(age_out_celltype[[5]]$gene%in%c("AGR2","EZR","LGALS2", "LGALS3", "ITGB1", "RER1")),]




### smooth line violin plot
data_plot<-FetchData(seu_cleaned_age, vars = c("AGR2","EZR","LGALS2", "LGALS3", "ITGB1", "RER1"))
data_plot$cell_id<-rownames(data_plot)
data_plot<-melt(data_plot)
meta<-seu_cleaned_age@meta.data
meta$cell_id<-rownames(meta)

data_plot<-merge(data_plot, meta, by="cell_id")
colnames(data_plot)[which(colnames(data_plot)=="variable")]<-"gene"

Means <- data_plot %>% group_by(age_con, gene,Gamma.Annotation) %>%
  summarize(Avg = mean(value))

age_cholangiocyte<-ggplot() +
  geom_violin(data = data_plot, mapping = aes(x = age_con, y = value), color="grey80", fill="grey80") +
  #geom_point(data = Means, mapping = aes(x = age_con, y = Avg), color="red") +
  geom_line(data = Means, mapping = aes(x = age_con, y = Avg, group=gene), color="black")+
  #geom_smooth(data = Means, mapping = aes(x = age_con, y = Avg, group=gene), method="lm")+
  facet_grid(gene~Gamma.Annotation)+theme_bw()+xlab("Age")+ylab("Expression")+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(data = sig_genes_relaxed[which(sig_genes_relaxed$gene%in%c("AGR2","EZR","LGALS2", "LGALS3", "ITGB1", "RER1")),],
            aes(x = 2, label = round(coef,2), y = ceiling(max(data_plot$value))), size=3)+
  ylim(0, ceiling(max(data_plot$value))*1.1)
age_cholangiocyte

ggsave(age_cholangiocyte, file="scratch/HLiCA/age_cholangiocyte_genes.png", h=8, w=8)


