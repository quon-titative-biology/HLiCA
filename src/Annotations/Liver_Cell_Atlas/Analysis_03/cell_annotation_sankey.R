library(ggplot2)
library(here)
library(dplyr)
library(scales)
library(reshape2)
library(cowplot)
library(ggsankey)


source("scripts/00_pretty_plots.R")
source("scripts/colour_palette.R")


beta_annotations<-read.csv("data/02_Beta_Annotations.csv")

load(here("data/meta_amalgamation.RData"))


meta_HLiCA$cell<-rownames(meta_HLiCA)

nrow(beta_annotations)
nrow(meta_HLiCA)

nrow(beta_annotations)-nrow(meta_HLiCA)

dim(meta_HLiCA[which(rownames(meta_HLiCA)%in%beta_annotations$X),])
dim(beta_annotations[which(!(beta_annotations$X%in%rownames(meta_HLiCA))),])

table(beta_annotations[which(!(beta_annotations$X%in%rownames(meta_HLiCA))),]$Potential.Doublets)
table(beta_annotations[which(!(beta_annotations$X%in%rownames(meta_HLiCA))),]$Beta.Annotation.Sublineage)

55431-43004
12427


table(beta_annotations[which(!(beta_annotations$X%in%rownames(meta_HLiCA))),]$Potential.Doublets.2)

missing_why<-beta_annotations[which(!(beta_annotations$X%in%rownames(meta_HLiCA))),]
missing_why<-missing_why[which(missing_why$Potential.Doublets==FALSE | is.na(missing_why$Potential.Doublets)),]
nrow(missing_why)
table(missing_why$Potential.Doublets)
table(missing_why$Alpha.Annotations)
table(missing_why$Beta.Annotation.SubLineage)

table(missing_why$Alpha.Annotations,missing_why$Beta.Annotation.SubLineage)
2408+846+5272
189+3712

#TBC 8526 


df_forsankey<-beta_annotations
df_forsankey$Alpha.Annotations.Doublets<-df_forsankey$Alpha.Annotations
df_forsankey$Alpha.Annotations.Doublets[which(df_forsankey$Potential.Doublets==TRUE)]<-"Potential Doublet"

df_forsankey$Beta.Annotations<-df_forsankey$Beta.Annotation.Lineage
df_forsankey$Beta.Annotations[which(df_forsankey$Beta.Annotation.SubLineage=="TBC")]<-"TBC"
df_forsankey$Beta.Annotations[which(df_forsankey$Potential.Doublets==TRUE)]<-"Potential Doublet"
df_forsankey$Beta.Annotations[which(df_forsankey$Alpha.Annotations=="Erythrocyte")]<-"Erythrocyte"
df_forsankey$Beta.Annotations[which(df_forsankey$Alpha.Annotations=="Doublet")]<-"Doublet"


meta_HLiCA$cell<-rownames(meta_HLiCA)
meta_HLiCA$lineage_alpha<-meta_HLiCA$Alpha.Annotations
df_forsankey<-merge(df_forsankey, meta_HLiCA[,c("cell","Gamma.Annotation","lineage_alpha")], by.x="X", by.y="cell", all.x = T)
table(df_forsankey$lineage_alpha, df_forsankey$Beta.Annotation.Lineage)
table(df_forsankey$Alpha.Annotations, df_forsankey$Beta.Annotation.Lineage)

df_forsankey$Gamma.Annotations<-df_forsankey$Alpha.Annotations
df_forsankey$Gamma.Annotations[which(df_forsankey$Potential.Doublets==TRUE)]<-"removed"
df_forsankey$Gamma.Annotations[which(is.na(df_forsankey$Gamma.Annotation))]<-"removed"
df_forsankey$Gamma.Annotations[which(df_forsankey$Gamma.Annotation=="pDC")]<-"pDC"
df_forsankey$Gamma.Annotations[which(df_forsankey$Gamma.Annotation%in%c("Cycling Cells","Cycling"))]<-"Cycling"
df_forsankey$Gamma.Annotations[which(df_forsankey$Gamma.Annotation=="doublets")]<-"removed"
df_forsankey$Gamma.Annotations[which(df_forsankey$Gamma.Annotation=="removed")]<-"removed"

save(df_forsankey, file=here("data/annotation_tracking.RData"))





load(here("data/annotation_tracking.RData"))
dim(df_forsankey)
head(df_forsankey)






#meta<-read.csv(here("/media/redgar/Seagate Portable Drive/HLiCA/annotation_review/healthy_metadata_gamma.csv"))
df_forsankey$Gamma.Annotation[which(df_forsankey$Beta.Annotation.SubLineage=="TBC")]<-"removed"
df_forsankey$Gamma.Annotation[which(df_forsankey$Alpha.Annotations=="Erythrocyte")]<-"removed"
df_forsankey$Gamma.Annotation[which(df_forsankey$Alpha.Annotations=="Doublet")]<-"removed"
df_forsankey$Beta.Annotation.SubLineage[which(df_forsankey$Alpha.Annotations=="Erythrocyte")]<-"removed"
df_forsankey$Beta.Annotation.SubLineage[which(df_forsankey$Alpha.Annotations=="Doublet")]<-"removed"
df_forsankey$Gamma.Annotation[which(df_forsankey$Beta.Annotation.SubLineage%in%c("Stressed Cells","Poor Quality","Doublet"))]<-"removed"
df_forsankey$Beta.Annotation.SubLineage[which(df_forsankey$Beta.Annotation.SubLineage%in%c("Doublet"))]<-"TBC"
df_forsankey$Gamma.Annotation[which(df_forsankey$Gamma.Annotation%in%c("doublets"))]<-"removed"



###############
## Bright dim fix and other name fixes
###############
df_forsankey$Gamma.Annotation[which(df_forsankey$Gamma.Annotation=="NRXN1+ Hepatic Stellate Cell")]<-"NRXN1+ Stromal Cell"
df_forsankey$Gamma.Annotation[which(df_forsankey$Gamma.Annotation=="Dim NK Cell")]<-"Dim NK Cell wrong"
df_forsankey$Gamma.Annotation[which(df_forsankey$Gamma.Annotation=="Bright NK Cell")]<-"Dim NK Cell"
df_forsankey$Gamma.Annotation[which(df_forsankey$Gamma.Annotation=="Dim NK Cell wrong")]<-"Bright NK Cell"
df_forsankey$Gamma.Annotation[which(df_forsankey$Gamma.Annotation=="LAMC2+ Small")]<-"LAMC2+ Cholangiocyte"
df_forsankey$Gamma.Annotation[which(df_forsankey$Gamma.Annotation=="Small ApoLipo")]<-"ApoLipo"
df_forsankey$Gamma.Annotation[which(df_forsankey$Gamma.Annotation=="Small Keratin")]<-"Keratin"
df_forsankey$Gamma.Annotation[which(df_forsankey$Gamma.Annotation=="CXCL8+ Small Keratin")]<-"CXCL8+ Keratin"
df_forsankey$Gamma.Annotation[which(df_forsankey$Gamma.Annotation=="Large Mucus Secreting")]<-"Mucus Secreting"




df <- df_forsankey %>%
  make_long(Alpha.Annotations,Beta.Annotation.SubLineage,Gamma.Annotation)


df_forsankey$Gamma.Annotation[which(df_forsankey$Beta.Annotation.SubLineage=="Potentially Technical")]
df_forsankey$Gamma.Annotation[which(df_forsankey$Beta.Annotation.SubLineage=="Vascular Smooth Muscle Cell")]
df_forsankey$Gamma.Annotation[which(df_forsankey$Gamma.Annotation=="doublets")]

celltype_order<- c( "Cycling","Cycling Cells","pDC", "pDCs",
                   "Myeloid Cell",
                   "Kupffer Cells", "Classical Monocytes", "Type 2 cDCs","Neutrophils",
                   "Non-classical Monocytes","Non Classical Monocytes", "TREM2+ Macrophages", "MAMLD1+ Trans Monocytes","LAM-Like",
                   "DCmac","Macrophages (cDCs?)",  "Activated Monocytes"," Activated Macrophages", "migDCs", "Type 1 cDCs", "Type 1 cDCs ",
                   "BAG3+ Monocytes", "T Cells?",
                   
                   "Cholangiocyte",
                   "Keratin / Small / Ribosomal / Cycling","ApoLipo / Small",  "LAMC2+ / Small","Small (Keratin)",     "ApoLipo",   
                   "Stressed?" ,"Keratin", "LAMC2+ Cholangiocyte","CXCL8+ Keratin",
                   "Mucus Secreting / Large / Keratin","Keratin / CXCL8+" , "Mucus Secreting", 
                   
                   "Endothelia",
                    "Portal Vein",  "Central Vein" ,"Vascular Endothelial Cell",
                   "Lymphatic", "Hepatic Artery","Periportal LSEC","Central Venous LSEC", "LSEC", 
                   
                   "Mesenchyme", 
                   "Portal Fibroblast", "CUX2+ Hepatic Stellate Cell", "NRXN1+ Stromal Cell","Hepatic Stellate Cell" ,"Vascular Smooth Muscle Cell",
                   
                   "Lymphocyte",
                   "Plasma B Cell","Naive B Cell", "B Cell", "IgA B cells", "IgG B cells",
                   
                   "Bright NK Cell", "cNK / Bright NK Cell"   , "Dim NK Cell", "Liver Resident / Dim NK Cell",
                   
                   "Helper T Cell","CD4+  / Helper T Cell" , "Exhausted CD8 T Cell","CD8 T Cell","CD8+ / Cytotoxic T Cell", 
                   "MAIT T Cell", "CD4+ / MAIT T Cell", "Regulatory T Cell","CD4+ / Regulatory T Cell"  , 
                   "Stressed Cells",  
                   "Poor Quality",
                   
                   "Hepatocyte",
                   "Ribosomal+ Hepatocyte","Ribosomal+ / APOA2+"  ,"UGT+ Hepatocyte" ,"Mito+ Hepatocyte" ,
                   "Pericentral Hepatocyte", "Periportal Hepatocyte","SERPINE1+ Hepatocyte",
                   "ORM2+ / APOE+","Ribosomal+ / FTL+" ,   
                   "Central","Portal",
                   "SERPINE1+","Mito+","Potentially Technical" ,
                   
                   "Erythrocyte","TBC","removed","Doublet"
)


unique(df_forsankey$Alpha.Annotations[which(!(df_forsankey$Alpha.Annotations%in%celltype_order))])
unique(df_forsankey$Beta.Annotation.SubLineage[which(!(df_forsankey$Beta.Annotation.SubLineage%in%celltype_order))])
unique(df_forsankey$Gamma.Annotation)[which(!(unique(df_forsankey$Gamma.Annotation)%in%celltype_order))]


df$node <- factor(df$node,levels = rev(celltype_order))

df$next_node <- factor(df$next_node,levels = rev(celltype_order))

levels(df$x)<-c("Initial annotation" ,"Annotation shared\nwith bionetwork\n for input", "Final Annotation") 

sankey_all_annotation<-ggplot(df, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey() +
  theme_sankey(base_size = 16)+  geom_sankey(flow.alpha = 0.75, node.color = 1,linewidth=0.01) +
  geom_sankey_label(aes(label=factor(node)),size = 2.25, color = 1, fill = "white",hjust = 0)+fillscale_sankey+
  theme(legend.position = "none")+
  xlab("")
save_plts(sankey_all_annotation, "sankey_all_annotation_alphabetagamma", h=10, w=16)


sankey_all_annotation<-ggplot(df, aes(x = x, 
                                      next_x = next_x, 
                                      node = node, 
                                      next_node = next_node,
                                      fill = factor(node))) +
  geom_sankey() +
  theme_sankey(base_size = 10)+  geom_sankey(flow.alpha = 0.75, node.color = 1,linewidth=0) +
  geom_sankey_label(aes(label=factor(node)),size = 1.25, color = 1.5, fill = "white",hjust = 0,label.size = 0.1)+fillscale_sankey+
  theme(legend.position = "none")+
  xlab("")
save_plts(sankey_all_annotation, "sankey_all_annotation_alphabetagamma_small", h=5.3, w=5.3)

