library(tiff)
library(here)

#remotes::install_github(repo = "satijalab/seurat", ref = "visium-hd")
#https://satijalab.org/seurat/articles/visiumhd_commands_intro
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(arrow)
library(reshape2)
library(cowplot)
library(ggforce)
library(colorspace)

source(here("scripts/00_pretty_plots.R"))

UGThepatocyte<-read.csv(here("../HLiCA/data/DEG/Hepatocytes_UGT+ Hepatocyte.csv"))
UGThepatocyte<-UGThepatocyte[which(UGThepatocyte$avg_log2FC>0),]






### load VisiumHD
annotation_C107_8um<-read.csv(here("data/metadata_for_C107.csv"))

localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

identical(annotation_C107_8um$X, colnames(object))
rownames(annotation_C107_8um) <- annotation_C107_8um$X
annotation_C107_8um$X<-NULL

object<-AddMetaData(object, annotation_C107_8um)

object$cell<-colnames(object)

# Prepare spot coordinates
coords <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))

# want to show gene expression
counts <- object@assays$Spatial$counts



rownames(counts)[grep("UGT", rownames(counts))]

UGThepatocyte<-read.csv(here("../HLiCA/data/DEG/Hepatocytes_UGT+ Hepatocyte.csv"))
UGThepatocyte<-UGThepatocyte[which(UGThepatocyte$avg_log2FC>0),]
rownames(counts)[which(rownames(counts)%in%UGThepatocyte$gene)]

DotPlot(object, features=c(rownames(counts)[grep("UGT", rownames(counts))], UGThepatocyte$gene), group.by = "final_anno")

### will use most highly expressed UGT, plus UGT cell type marker: CYP26A1 and sex differential gene: CYP3A4

# Flip y for coords too
# coords$y_centred <- (coords$y - ymin)+1
# coords$x_centred <- (xmin - (coords$x)+1)+(xmax-xmin)

bin_size <- 13.4543/4
half_bin <- bin_size / 2
unit_per_um <- bin_size / 2

# Choose scale bar size in μm
scale_length_um <- 1000
scale_length_plot <- scale_length_um * unit_per_um
x_start <- (-max(coords$x))
y_start <- max(coords$y)









########################
## Gene expression plot
########################
# UGT2B4 UGT2B10 UGT1A1 CYP3A4
genes<-c("CYP26A1","UGT2B4", "UGT2B10", "UGT1A1","CYP3A4") 
gene_exp<-FetchData(object, vars= genes)
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)
annotation_C107_8um$cell<-rownames(annotation_C107_8um)
plt<-merge(gene_exp, annotation_C107_8um, by="cell")

plt<-merge(plt, coords[,c("cell","x","y")], by="cell")



## summarize
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_hep<-plt[grep(c("hepatocyte"), plt$final_anno),]

plt_summary<-plt_hep %>% group_by(final_anno, variable) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
plt_summary <- plt_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(mn))
plt_summary<-as.data.frame(plt_summary)


# remove dots where 0 cell expressing marker
plt_summary<-plt_summary[(which(plt_summary$count>0)),]

fancy_dotplot<-plot_grid(
  ggplot(plt_summary, aes(final_anno, variable, color=scaled, size=percent_exp))+geom_point()+
    theme_classic()+
    scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
    scale_size(name="Percent\nCells\nExpressing")+
    theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank(),
          plot.margin = margin(0.25,0.25,0,0.25,"cm"))+
    ggplot(plt_summary, aes(final_anno, y=1, fill=final_anno))+geom_tile(color="black")+
    theme_classic()+fillscale_cellType+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",axis.line  = element_blank(),
          plot.margin = margin(0,0,1,1,"cm")),
  ncol=1, rel_heights = c(1.3,1), align = "v", axis="lr")
fancy_dotplot


table(plt_hep$value, plt_hep$variable)

# these are very lowly expressed genes
# 99% 95% 99% 93% of cell expressing the transcript only have a count of 1
# so makes more sense to show the change in %bin expressing I think



plt_summary$final_anno<-as.factor(plt_summary$final_anno)
levels(plt_summary$final_anno)<-1:9

# layer_expression<-ggplot(plt_summary, aes(final_anno, percent_exp, color=final_anno, fill=final_anno, group=variable))+
#   geom_line(size=1)+
#   geom_point(aes(size=scaled), shape=21, color="white")+
#   geom_text(aes(label=variable, hjust = -0.3),plt_summary[which(plt_summary$final_anno==1),], color="black")+
#   theme_bw()+ylab("Percent Cells Expressing")+xlab("Lobule Layer")+
#   scale_color_manual(values = rev(c("#184E77", "#1E6091","#1A759F", "#168AAD", "#34A0A4", "#52B69A", "#76C893", "#99D98C" ,"#B5E48C")))+
#   scale_fill_manual(values = rev(c("#184E77", "#1E6091","#1A759F", "#168AAD", "#34A0A4", "#52B69A", "#76C893", "#99D98C" ,"#B5E48C")))+
#   scale_size_continuous(name = "Scaled\nMean\nExpression") 
# layer_expression
# save_plts(layer_expression, "layer_expression_UGT_forHLiCA", w=4,h=6)                         



layer_expression<-ggplot(plt_summary, aes(final_anno, mn, color=final_anno, fill=final_anno, group=variable))+
  geom_line(size=1)+
  geom_point(aes(size=percent_exp), shape=21, color="white")+
  geom_text(aes(label=variable, hjust = -0.3),plt_summary[which(plt_summary$final_anno==1),], color="black")+
  theme_bw()+ylab("Mean Expression")+xlab("Lobule Layer")+
  scale_color_manual(values = rev(c("#184E77", "#1E6091","#1A759F", "#168AAD", "#34A0A4", "#52B69A", "#76C893", "#99D98C" ,"#B5E48C")))+
  scale_fill_manual(values = rev(c("#184E77", "#1E6091","#1A759F", "#168AAD", "#34A0A4", "#52B69A", "#76C893", "#99D98C" ,"#B5E48C")))+
  scale_size_continuous(name = "Percent\nCells\nExpressing ") 
layer_expression
save_plts(layer_expression, "layer_expression_UGT_forHLiCA", w=4,h=6)                         







merscope_UGT<-plot_grid(
  ggplot(MERSCOPE_UGT[rev(order(MERSCOPE_UGT$Cell_Type)),], aes(x,y, color=Cell_Type))+
    geom_point(size=0.2)+facet_wrap(~sample_id, scales = "free")+theme_bw()+
    scale_color_manual(values = c("#FFD700","#8B4513","#4b911d","forestgreen","#08519c","#6495ED","#9ecae1","#b80783","#914e84"))+ 
    guides(color = guide_legend(override.aes = list(size = 2))),
  ggplot(plt[which(plt$Cell_Type%in%c("Hep_1","Hep_2","Hep_3")),], aes(Cell_Type,log(value)))+
    geom_violin(fill="grey90",color="grey90")+geom_boxplot(aes(fill=Cell_Type),width=0.5,outlier.shape=NA)+
    facet_wrap(~variable)+theme_bw()+
    scale_fill_manual(values = c("#08519c","#6495ED","#9ecae1"))+
    theme(legend.position = "none"), rel_widths = c(3,1))
save_plts(merscope_UGT, "UGT_MERFISH", w=14, h=4)




#################
## hepaocyte zonation plot
#################
bin_size <- 13.4543
half_bin <- bin_size / 2
unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 1000
scale_length_plot <- scale_length_um * unit_per_um
x_start <- (-max(plt$x))
y_start <- max(plt$y)

get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

plt$final_anno<-as.factor(plt$final_anno)
levels(plt$final_anno)<-gsub("bin","layer", levels(plt$final_anno))

plot_bins<-ggplot() +
  geom_point(data = plt,aes(x = y, y = -x), color = "black",size = 0.8) +
  geom_point(data = plt,aes(x = y, y = -x), color = "grey90",size = 0.1) +
  geom_point(data = plt[grep("hepatocyte",plt$final_anno),],aes(x = y, y = -x, color = final_anno),size = 0.1, alpha=0.1) +
  annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start+ 300,yend = x_start+300,color = "black",linewidth = 1) +
  annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start + 200,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
  coord_fixed() +colscale_cellType+
  theme_void()
plot_bins

full_tissue<-plot_grid(plot_bins+theme(legend.position = "none"), get_leg(plot_bins+guides(color = guide_legend(ncol=4,override.aes = list(size = 2)))), ncol=1, rel_heights = c(3,1))
full_tissue

save_plts(full_tissue, "HLiCA_C107_fulltissue_tissehepatocytes_only", w=8, h=8)





