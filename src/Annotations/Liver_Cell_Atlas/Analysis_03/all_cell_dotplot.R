### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(cowplot)
library(RColorBrewer)
library(colorspace)


source("scripts/00_pretty_plots.R")
source("scripts/colour_palette.R")



####################################
## Fancy dot plot
####################################
d10x<-readRDS(file="../healthy_RNA_merged_harmonized_gamma.rds")

DefaultAssay(d10x)<-"RNA"

## markers
myeloid<-c("CD68","CD163","MARCO","C1QC","TREM2","LAMP3","CD9","PLXDC2","BAG3","MAMLD1","CLEC9A","CLEC10A","LYZ","IL1B","FCGR3A","FCGR3B")

lymphocyte<-c("CD8A","CD8B","SLC4A10","KLRF1","FCER1G","IL7R","FOXP3","PDCD1","MS4A1","CD19","JCHAIN","IGHA1","IGHG1","GZMB")

cholangiocyte<-c("EPCAM","KRT18","KRT19","APOA2","APOA1","MUC5B","LAMC2","CXCL8")
mesenchyme<-c("COL1A1","SPON1","RGS5","RGS6","RELN", "CUX2","NRXN1")
endothelial<-c("PECAM1","CD34","CPE","TAGLN","CD36","LYVE1","CLEC4G","PDPN")

hepatocytes<-c("ALB", "HAL", "CYP3A4", "CYP2E1","RPS6","UGT2B10","UGT2B7","MT−RNR2","MT−ND3","SERPINE1")

cycling<-c("MKI67","TOP2A")



gene_exp<-FetchData(d10x, vars=c(myeloid,lymphocyte,cholangiocyte,mesenchyme,endothelial,hepatocytes,cycling))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)
d10x@meta.data$cell<-rownames(d10x@meta.data)

plt<-merge(gene_exp, d10x@meta.data, by="cell")

## summarize
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_summary<-plt %>% group_by(Gamma.Annotation, variable) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
plt_summary <- plt_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(mn))
plt_summary<-as.data.frame(plt_summary)

# remove dots where 0 cell expressing marker
plt_summary<-plt_summary[(which(plt_summary$count>0)),]
plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(myeloid,lymphocyte,cholangiocyte,mesenchyme,endothelial,hepatocytes,cycling)))

###############
## Bright dim fix and other name fixes
###############
plt_summary$Gamma.Annotation[which(plt_summary$Gamma.Annotation=="NRXN1+ Hepatic Stellate Cell")]<-"NRXN1+ Stromal Cell"
plt_summary$Gamma.Annotation[which(plt_summary$Gamma.Annotation=="Dim NK Cell")]<-"Dim NK Cell wrong"
plt_summary$Gamma.Annotation[which(plt_summary$Gamma.Annotation=="Bright NK Cell")]<-"Dim NK Cell"
plt_summary$Gamma.Annotation[which(plt_summary$Gamma.Annotation=="Dim NK Cell wrong")]<-"Bright NK Cell"
plt_summary$Gamma.Annotation[which(plt_summary$Gamma.Annotation=="LAMC2+ Small")]<-"LAMC2+ Cholangiocyte"
plt_summary$Gamma.Annotation[which(plt_summary$Gamma.Annotation=="Small ApoLipo")]<-"ApoLipo"
plt_summary$Gamma.Annotation[which(plt_summary$Gamma.Annotation=="Small Keratin")]<-"Keratin"
plt_summary$Gamma.Annotation[which(plt_summary$Gamma.Annotation=="CXCL8+ Small Keratin")]<-"CXCL8+ Keratin"
plt_summary$Gamma.Annotation[which(plt_summary$Gamma.Annotation=="Large Mucus Secreting")]<-"Mucus Secreting"



plt_summary$Gamma.Annotation<-factor(plt_summary$Gamma.Annotation, levels=(c(  "Kupffer Cells", "TREM2+ Macrophages","migDCs", "DCmac",
                                                                               "BAG3+ Monocytes", "MAMLD1+ Trans Monocytes", 
                                                                               "Type 1 cDCs","Type 2 cDCs","Classical Monocytes", "Activated Monocytes", 
                                                                               "Non-classical Monocytes",  "Neutrophils", 
                                                                               
                                                                               "CD8 T Cell", "MAIT T Cell", "Bright NK Cell", "Dim NK Cell", "Helper T Cell", 
                                                                               "Regulatory T Cell", "Exhausted CD8 T Cell", 
                                                                               "B Cell", "Plasma B Cell","IgA B cells", "IgG B cells","pDC",
                                                                               
                                                                               "Keratin", "ApoLipo", "Mucus Secreting", 
                                                                               "LAMC2+ Cholangiocyte","CXCL8+ Keratin",
                                                                               
                                                                               "Portal Fibroblast","Vascular Smooth Muscle Cell",
                                                                               "Hepatic Stellate Cell","CUX2+ Hepatic Stellate Cell", 
                                                                               "NRXN1+ Stromal Cell",
                                                                               
                                                                               "Hepatic Artery","Portal Vein","Vascular Endothelial Cell", 
                                                                               "Periportal LSEC","Central Venous LSEC",  "LSEC", "Lymphatic", 
                                                                               
                                                                               "Periportal Hepatocyte","Pericentral Hepatocyte", "Ribosomal+ Hepatocyte",
                                                                               "UGT+ Hepatocyte" ,"Mito+ Hepatocyte" , "SERPINE1+ Hepatocyte",
                                                                               
                                                                               "Cycling")))

fancy_dotplot<-plot_grid(
  ggplot(plt_summary, aes(Gamma.Annotation, variable, color=scaled, size=percent_exp))+geom_point()+
    th+theme_classic()+
    scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
    scale_size(name="Percent\nCells\nExpressing")+
    theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank(),
          plot.margin = margin(0.25,0.25,0,0.25,"cm"))+
    geom_hline(yintercept = c(2.5, 10.5, 18.5,25.5,33.5,34.5,47.5), color="grey70")+
    geom_vline(xintercept = c(12.5,23.5,24.5,29.5,34.5,41.5), color="grey70"),
  ggplot(plt_summary, aes(Gamma.Annotation, y=1, fill=Gamma.Annotation))+geom_tile(color="black")+
    th+theme_classic()+fillscale_cellType+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",axis.line  = element_blank(),
          plot.margin = margin(0,0,1,1,"cm")),
  ncol=1, rel_heights = c(4,0.7), align = "v", axis="lr")
fancy_dotplot
ggsave(file="figures/HLiCA_dot_plot_celltype.pdf", w=12,h=15)
