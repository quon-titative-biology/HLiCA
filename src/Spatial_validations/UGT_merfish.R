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
library(SCINA)
library(sceasy)

source(here("scripts/00_pretty_plots.R"))


schwann_markers<-read.csv(here("data/DEG/NRXN1+ Hepatic Stellate Cell.csv"))
schwann_markers<-schwann_markers[which(schwann_markers$avg_log2FC>0),]

novel_cell_markers<-c("NRXN1","MAMLD1","TREM2","BAG3",
                      "MUC5B","MUC13","MUC1","MUC3A","MUC12","MUC4","MUC20","TCN1","FCGBP","ARG1",
                      "SERPINE1","SGMS2","UGT2B7", "UGT1A7","CYP2A7P1")

MUC_markers<-c("MUC5B","MUC13","MUC1","MUC3A","MUC12","MUC4","MUC20","TCN1","FCGBP","ARG1")

UGThepatocyte<-read.csv(here("data/DEG/Hepatocytes_UGT+ Hepatocyte.csv"))
UGThepatocyte<-UGThepatocyte[which(UGThepatocyte$avg_log2FC>0),]

################
## MERFISH
################
#https://datadryad.org/dataset/doi:10.5061/dryad.37pvmcvsg
MERFISH <- convertFormat(here("/media/redgar/Seagate Portable Drive/HLiCA/public_spatial/adata_healthy_merfish.h5ad"), from="anndata", to="seurat", main_layer="counts")
novel_cell_markers[which(novel_cell_markers%in%rownames(MERFISH))]
schwann_markers$gene[which(schwann_markers$gene%in%rownames(MERFISH))]
UGThepatocyte$gene[which(UGThepatocyte$gene%in%rownames(MERFISH))]

rownames(MERFISH)[grep("UGT",rownames(MERFISH))]
rownames(MERFISH)[grep("CYP",rownames(MERFISH))]

################
## UGT (MERFISH)
################
DimPlot(MERFISH, group.by = "Cell_Type")

DimPlot(MERFISH, group.by = "Cell_Type", reduction = "spatial", split.by = "sample_id")+
  scale_color_manual(values = c("#FF8C00","#8B4513","#FFD700","#F0E68C","#08519c","#6495ED","#9ecae1","#b80783","#914e84"))
FeaturePlot(MERFISH, features = "SERPINE1", reduction = "spatial", split.by = "sample_id")
FeaturePlot(MERFISH, features = "SGMS2", reduction = "spatial", split.by = "sample_id")

plot_grid(DimPlot(MERFISH, group.by = "Cell_Type", reduction = "spatial", split.by = "sample_id")+
            scale_color_manual(values = c("#FF8C00","#8B4513","#FFD700","#F0E68C","#08519c","#6495ED","#9ecae1","#b80783","#914e84")),
          FeaturePlot(MERFISH, features = "SERPINE1", reduction = "spatial", split.by = "sample_id"), ncol=1)


DotPlot(MERFISH, features = c("SERPINE1","SGMS2"), group.by = "Cell_Type")

FeaturePlot(MERFISH, features = c("UGT2B10","UGT2B4"))


MERFISH_UGT<-MERFISH@meta.data

gene_exp<-FetchData(MERFISH, vars= c("UGT2B10","UGT2B4","CYP3A4"))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)
MERFISH_UGT$cell<-rownames(MERFISH_UGT)
plt<-merge(gene_exp, MERFISH_UGT, by="cell")


plt_hep<-plt[which(plt$Cell_Type%in%c("Hep_1","Hep_2","Hep_3")),]
plt_hep$Cell_Type<-factor(plt_hep$Cell_Type, levels=c("Hep_3","Hep_2","Hep_1"))

MERFISH_UGT_plt<-ggplot()+
  geom_point(data=MERFISH_UGT, aes(x,y), color="black",size=1)+
  geom_point(data=MERFISH_UGT, aes(x,y), color="grey90", size=0.2)+
  geom_point(data=MERFISH_UGT[grep("Hep", MERFISH_UGT$Cell_Type),], aes(x,y, color=Cell_Type), size=0.2)+
  theme_void()+
  scale_color_manual(values = c("#184E77","#34A0A4","#99D98C"))+ 
  facet_wrap(~sample_id, scales = "free")+guides(color = guide_legend(override.aes = list(size = 2)))
MERFISH_UGT_plt
save_plts(MERFISH_UGT_plt, "UGT_MERFISH", w=12, h=3)


## plot one sample only
merfish_AM042<-MERFISH_UGT[which(MERFISH_UGT$sample_id=="AM042"),]

merfish_AM042_plt<-ggplot()+
  geom_point(data=merfish_AM042, aes(x,y), color="black",size=1)+
  geom_point(data=merfish_AM042, aes(x,y), color="grey90", size=0.2)+
  geom_point(data=merfish_AM042[grep("Hep", merfish_AM042$Cell_Type),], aes(x,y, color=Cell_Type), size=0.2)+
  theme_bw()+
  scale_color_manual(values = c("#184E77","#34A0A4","#99D98C"))+ 
  guides(color = guide_legend(override.aes = list(size = 2)))+theme_void()
merfish_AM042_plt
save_plts(merfish_AM042_plt+theme(legend.position = "none"), "UGT_MERFISH_AM042", w=4, h=4)

get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

save_plts(get_leg(merfish_AM042_plt), "UGT_MERFISH_AM042_leg", w=3, h=3)

UGT_violin<-ggplot(plt_hep, aes(Cell_Type,log(value)))+
  geom_violin(fill="grey90",color="grey90")+geom_boxplot(aes(fill=Cell_Type),width=0.5,outlier.shape=NA)+
  facet_wrap(~variable)+theme_bw()+
  scale_fill_manual(values = c("#9ecae1","#6495ED","#08519c"))+
  theme(legend.position = "none")
UGT_violin
save_plts(UGT_violin, "UGT_expression_merfish", w=4, h=4)


### DEG
Idents(MERFISH)<-"Cell_Type"

deg13<-FindMarkers(MERFISH, ident.1 = "Hep_1", ident.2 = "Hep_3",logfc.threshold = 0)
deg13[grep("UGT", rownames(deg13)),]

deg23<-FindMarkers(MERFISH, ident.1 = "Hep_2", ident.2 = "Hep_3",logfc.threshold = 0)
deg23[grep("UGT", rownames(deg23)),]



#### plot as perfect expressing
plt_summary<-plt_hep %>% group_by(Cell_Type, variable) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)





# UGT are these are very lowly expressed genes
ggplot(plt_summary, aes(Cell_Type, percent_exp, color=Cell_Type, fill=Cell_Type, group=variable))+
  geom_line(size=1)+
  geom_point(aes(size=mn), shape=21, color="white")+
  geom_text(aes(label=variable, hjust = -0.3),plt_summary[which(plt_summary$Cell_Type=="Hep_3"),], color="black")+
  theme_bw()+ylab("Percent Cells Expressing")+xlab("Lobule Layer")+
  scale_fill_manual(values = c("#184E77","#34A0A4","#99D98C"))+ 
  scale_color_manual(values = c("#184E77","#34A0A4","#99D98C"))+ 
  scale_size_continuous(name = "Scaled\nMean\nExpression") 

# scale_color_manual(values = c("#9ecae1","#6495ED","#08519c"))+
#   scale_fill_manual(values = c("#9ecae1","#6495ED","#08519c"))+

layer_expression<-ggplot(plt_summary, aes(Cell_Type, mn, color=Cell_Type, fill=Cell_Type, group=variable))+
  geom_line(size=1)+
  geom_point(aes(size=percent_exp), shape=21, color="white")+
  geom_text(aes(label=variable, hjust = -0.3),plt_summary[which(plt_summary$Cell_Type=="Hep_3"),], color="black")+
  theme_bw()+ylab("Mean Expression")+xlab("Lobule Layer")+
  scale_fill_manual(values = rev(c("#184E77","#34A0A4","#99D98C")))+ 
  scale_color_manual(values =rev(c("#184E77","#34A0A4","#99D98C")))+ 
  scale_size_continuous(name = "Percent\nCells\nExpressing") 
layer_expression
save_plts(layer_expression, "layer_expression_UGT_forHLiCA_merfish", w=4,h=6)                         
