library(ggplot2)
library(scales)
library(dplyr)
library(viridis)
library(cowplot)
library(here)
library(Seurat)
library(tiff)
library(reshape2)
library(ggforce)
library(colorspace)

source(here("scripts/00_pretty_plots.R"))


annotation_C107_8um<-read.csv(here("data/metadata_for_C107.csv"))

localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))



MUC_markers<-c("MUC5B","MUC13","MUC1","MUC3A","MUC12","MUC4","MUC20","TCN1","FCGBP","AGR2")


cholangiocyte_bins<-annotation_C107_8um[which(annotation_C107_8um$final_anno=="intrahepatic cholangiocyte"),]
cholangiocyte_object <- subset(object, cells = cholangiocyte_bins$X)
#rm(object)

cholangiocyte_object <- NormalizeData(cholangiocyte_object)
cholangiocyte_object <- FindVariableFeatures(cholangiocyte_object)
cholangiocyte_object <- ScaleData(cholangiocyte_object)
cholangiocyte_object <- RunPCA(cholangiocyte_object, reduction.name = "pca.008um")
cholangiocyte_object <- FindNeighbors(cholangiocyte_object, reduction = "pca.008um", dims = 1:30)
cholangiocyte_object <- FindClusters(cholangiocyte_object, resolution = 1, cluster.name = "seurat_cluster.008um")
cholangiocyte_object <- RunUMAP(cholangiocyte_object,reduction = "pca.008um",  reduction.name = "umap.008um", dims = 1:30)


DimPlot(cholangiocyte_object)
DotPlot(cholangiocyte_object, features = MUC_markers)
FeaturePlot(cholangiocyte_object, feature=c("AGR2","MUC3A","MUC5B","SERPINA1"))
FeaturePlot(cholangiocyte_object, feature=MUC_markers[1:4])
FeaturePlot(cholangiocyte_object, feature=MUC_markers[5:8])
FeaturePlot(cholangiocyte_object, feature=c("DSCAM","SEMA6D","KRT19","EPCAM"))

FeaturePlot(cholangiocyte_object, feature=c("AGR2","MUC3A","MUC5B","SERPINA1","MUC20","MUC13","MUC1","EPCAM"))
DotPlot(cholangiocyte_object, feature=c("MUC3A","MUC5B","MUC20","MUC13","MUC1","EPCAM","DSCAM","SEMA6D"))
DotPlot(cholangiocyte_object, feature=c("PHGR1","TFF1","TFF2","LYZ","MUC1","EPCAM","DSCAM","SEMA6D"))


de_0<-FindMarkers(cholangiocyte_object, ident.1="0")
de_0[grep("AGR|MUC", rownames(de_0)),]


MUC_cholangiocyte<-colnames(cholangiocyte_object)[which(cholangiocyte_object$seurat_clusters==0)]
save(MUC_cholangiocyte, file=here("data/mucus_secreting.RData"))

### Fancy dotplot
gene_exp<-FetchData(object, vars=MUC_markers)
gene_exp$cell<-rownames(gene_exp)
gc()

plt<-merge(gene_exp, annotation_C107_8um, by.x="cell", by.y="X")
gc()



############ plot gene exp
# Prepare spot coordinates
coords <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords$bin<-rownames(coords)
rm(object)
gc()

plt_annotation_spatial<-merge(coords,plt, by.x="bin", by.y="cell")

bin_size <- 13.4543
half_bin <- bin_size / 2
unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 1000
scale_length_plot <- scale_length_um * unit_per_um
x_start <- (-max(plt_annotation_spatial$x))
y_start <- max(plt_annotation_spatial$y)

get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

# plot_bins<-ggplot() +
#   geom_point(data = plt_annotation_spatial[which(plt_annotation_spatial$ARG1==0),],aes(x = y, y = -x),size = 0.25,color = "grey90") +
#   geom_point(data = plt_annotation_spatial[which(plt_annotation_spatial$ARG1!=0),],aes(x = y, y = -x, color = ARG1),size = 0.25) +
#   annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start+ 300,yend = x_start+300,color = "black",linewidth = 1) +
#   annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start + 200,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
#   geom_rect(aes(xmin = min(plt_annotation_spatial$y), xmax = max(plt_annotation_spatial$y), ymin = (-min(plt_annotation_spatial$x)), ymax = (-max(plt_annotation_spatial$x))),
#             fill = NA,color = "black",linewidth = 1) +
#   scale_color_viridis_c()+
#   theme_void()
# 
# 
# full_tissue<-plot_grid(plot_bins+theme(legend.position = "none"), get_leg(plot_bins+guides(color = guide_legend(ncol=4,override.aes = list(size = 2)))), ncol=1, rel_heights = c(3,1))
# full_tissue

plt_annotation_spatial$muc<-NA
plt_annotation_spatial$muc[which(plt_annotation_spatial$bin%in%MUC_cholangiocyte)]<-"Mucus Secreating"
table(plt_annotation_spatial$muc)

bin_size <- 13.4543
half_bin <- bin_size / 2
unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 1000
scale_length_plot <- scale_length_um * unit_per_um
x_start <- (-max(plt_annotation_spatial$x))
y_start <- max(plt_annotation_spatial$y)

get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

plt_annotation_spatial$final_anno<-as.factor(plt_annotation_spatial$final_anno)
levels(plt_annotation_spatial$final_anno)<-gsub("bin","layer", levels(plt_annotation_spatial$final_anno))


plot_bins<-ggplot() +
  geom_point(data = plt_annotation_spatial[which(is.na(plt_annotation_spatial$muc)),],aes(x = y, y = -x, color = final_anno),size = 0.1) +
  geom_rect(aes(xmin = min(plt_annotation_spatial$y), xmax = max(plt_annotation_spatial$y), ymin = (-min(plt_annotation_spatial$x)), ymax = (-max(plt_annotation_spatial$x))), fill = "white", alpha=0.5,linewidth = 1) +  
  geom_point(data = plt_annotation_spatial[which(!is.na(plt_annotation_spatial$muc)),],aes(x = y, y = -x),color="black",size = 0.1) +
  annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start+ 300,yend = x_start+300,color = "black",linewidth = 1) +
  annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start + 200,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
  geom_rect(aes(xmin = min(plt_annotation_spatial$y), xmax = max(plt_annotation_spatial$y), ymin = (-min(plt_annotation_spatial$x)), ymax = (-max(plt_annotation_spatial$x))), fill = NA,color = "black",linewidth = 1) +
  coord_fixed() +colscale_cellType+
  theme_void()
plot_bins


full_tissue<-plot_grid(plot_bins+theme(legend.position = "none"), get_leg(plot_bins+guides(color = guide_legend(ncol=4,override.aes = list(size = 2)))), ncol=1, rel_heights = c(3,1))
full_tissue

save_plts(full_tissue, "HLiCA_C107_fulltissue_tisse_MUC", w=8, h=8)





#plt_annotation_spatial_mini<-plt_annotation_spatial[sample(1:nrow(plt_annotation_spatial),100000),]

plt_annotation_spatial$muc<-NA
plt_annotation_spatial$muc[which(plt_annotation_spatial$bin%in%MUC_cholangiocyte)]<-"Mucus Secreating"
table(plt_annotation_spatial$muc)

plt_annotation_spatial$final_anno<-as.factor(plt_annotation_spatial$final_anno)
levels(plt_annotation_spatial$final_anno)<-gsub("bin","layer", levels(plt_annotation_spatial$final_anno))

xmin<-10200
xmax<-11500
ymin<-10200
ymax<-11600

plot_bins<-ggplot() +
  geom_point(data = plt_annotation_spatial,aes(x = y, y = -x), color = "black",size = 0.8) +
  geom_point(data = plt_annotation_spatial,aes(x = y, y = -x), color = "grey90",size = 0.1) +
  geom_point(data = plt_annotation_spatial[grep("hepatocyte",plt_annotation_spatial$final_anno),],
             aes(x = y, y = -x, color = final_anno),size = 0.1) +
  geom_point(data = plt_annotation_spatial[which(plt_annotation_spatial$final_anno%in%c("intrahepatic cholangiocyte")),],
             aes(x = y, y = -x, color = final_anno),size = 0.1) +
  geom_point(data = plt_annotation_spatial[which(!is.na(plt_annotation_spatial$muc)),],
             aes(x = y, y = -x),color="#eb150c",size = 0.1) +#FF00A0
  annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start+ 300,yend = x_start+300,color = "black",linewidth = 1) +
  annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start + 200,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
  coord_fixed() +colscale_cellType+
  geom_rect( aes(xmin = ymin, xmax = ymax, ymin = -xmin, ymax = -xmax),    fill = NA, color = "black", linewidth = 0.7)  +theme_void() 
plot_bins

full_tissue<-plot_grid(plot_bins+theme(legend.position = "none"), 
                       get_leg(plot_bins+guides(color = guide_legend(ncol=4,override.aes = list(size = 2)))), ncol=1, rel_heights = c(3,1))
full_tissue

save_plts(full_tissue, "HLiCA_C107_fulltissue_tisse_MUC_alternate", w=8, h=8)






# 
# 
# ggplot() +  geom_violin(data = plt_annotation_spatial,aes(y = final_anno, x =log10(SERPINE1), color = SERPINE1))
# 
# 
# hep_non0<-plt_annotation_spatial[which(plt_annotation_spatial$SERPINE1!=0),]
# 
# ggplot() +  geom_violin(data = hep_non0, aes(x = final_anno, y =(SERPINE1)))
# 
# table(plt_annotation_spatial$zones, plt_annotation_spatial$final_anno)
# 
# cholangiocyte_only<-plt_annotation_spatial[grep("intrahepatic cholangiocyte",plt_annotation_spatial$final_anno),]
# cholangiocyte_only$zones<-as.factor(cholangiocyte_only$zones)
# cholangiocyte_only$zones<-as.numeric(cholangiocyte_only$zones)
# ggplot(data = cholangiocyte_only, aes(x = zones, y =ARG1)) +  geom_point()+stat_smooth(method = "lm")
# 
# 
# 
# 
# ### dotplot
# gene_exp<-melt(plt_annotation_spatial[,c(1,32,which(colnames(plt_annotation_spatial)%in%c(novel_cell_markers,key_schwann)))])
# 
# 
# 
# ## summarize
# scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}
# 
# plt_summary<-gene_exp %>% group_by(final_anno, variable) %>% 
#   summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
# plt_summary <- plt_summary %>% group_by(variable) %>%
#   dplyr::mutate(scaled = scale_this(mn))
# plt_summary<-as.data.frame(plt_summary)
# 
# # remove dots where 0 cell expressing marker
# plt_summary<-plt_summary[(which(plt_summary$count>0)),]
# 
# plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(novel_cell_markers[2:19],key_schwann)))
# 
# plt_summary$final_anno<-factor(plt_summary$final_anno, levels=c(
#   "bin 1 hepatocyte", "bin 2 hepatocyte", "bin 3 hepatocyte", "bin 4 hepatocyte", "bin 5 hepatocyte", "bin 6 hepatocyte", "bin 7 hepatocyte",
#   "bin 8 hepatocyte", "bin 9 hepatocyte",
#   "endothelial cell of pericentral hepatic sinusoid","endothelial cell of periportal hepatic sinusoid",
#   "vein endothelial cell","endothelial cell of artery",
#   "fibroblast","hepatic stellate cell",
#   "intrahepatic cholangiocyte",
#   "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell",
#   "T cell","natural killer cell","hepatic pit cell",
#   "mature B cell","plasma cell", 
#   "neutrophil", "neutrophil2",
#   "monocyte","macrophage","Kupffer cell","mast cell","conventional dendritic cell", "plasmacytoid dendritic cell",
#   "erythrocyte", "Low_UMI","Unannotated","unknown"))
# 
# 
# fancy_dotplot<-plot_grid(
#   ggplot(plt_summary, aes(final_anno, variable, color=scaled, size=percent_exp))+geom_point()+
#     theme_classic()+
#     scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
#     scale_size(name="Percent\nCells\nExpressing")+
#     theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank())+
#     geom_hline(yintercept =c(2,3,10,12,14,18,19,20,22,25)+0.5, color="grey70")+
#     geom_vline(xintercept = c(9,13,15,16,21,23,25,31,32)+0.5, color="grey70"),
#   ggplot(plt_summary, aes(final_anno, y=1, fill=final_anno))+geom_tile(color="black")+
#     theme_classic()+fillscale_cellType+
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#           axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
#           legend.position = "none",axis.line  = element_blank(),
#           plot.margin = margin(t = 0,  # Top margin
#                                r = 50,  # Right margin
#                                b = 40,  # Bottom margin
#                                l = 10)),
#   ncol=1, rel_heights = c(6,2.5), align = "v", axis="lr")
# fancy_dotplot
