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
annotation_C107_8um$final_anno<-gsub("bin ","layer ",annotation_C107_8um$final_anno)


localdir <- "/media/redgar/Seagate Portable Drive/visiumHD_liver/MacParland_Diana__C107_D1/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

identical(annotation_C107_8um$X, colnames(object))
rownames(annotation_C107_8um) <- annotation_C107_8um$X
annotation_C107_8um$X<-NULL

object<-AddMetaData(object, annotation_C107_8um)

schwann_markers<-read.csv(here("../HLiCA/data/DEG/NRXN1+ Hepatic Stellate Cell.csv"))
schwann_markers<-schwann_markers[which(schwann_markers$avg_log2FC>0),]



DotPlot(object, features = c("NRXN1","PMP22","PLP1","COL1A1"), group.by = "final_anno")
DotPlot(object, features = c("NRXN1","PMP22","PLP1","COL1A1"), group.by = "spot_class")

endo_fibo<-rownames(annotation_C107_8um)[grep("endothe|fibro", annotation_C107_8um$final_anno)]

endo_fibo_object <- subset(object, cells = endo_fibo)

endo_fibo_object <- NormalizeData(endo_fibo_object)
endo_fibo_object <- FindVariableFeatures(endo_fibo_object)
endo_fibo_object <- ScaleData(endo_fibo_object)
endo_fibo_object <- RunPCA(endo_fibo_object, reduction.name = "pca.008um")
endo_fibo_object <- FindNeighbors(endo_fibo_object, reduction = "pca.008um", dims = 1:30)
endo_fibo_object <- FindClusters(endo_fibo_object, resolution = 0.6, cluster.name = "seurat_cluster.008um")
endo_fibo_object <- RunUMAP(endo_fibo_object,reduction = "pca.008um",  reduction.name = "umap.008um", dims = 1:30)


DimPlot(endo_fibo_object, label=T)
DotPlot(endo_fibo_object, features = c("NRXN1","PMP22","PLP1","COL1A1"))

FeaturePlot(endo_fibo_object, feature=c("NRXN1","PMP22","PLP1","COL1A1"))
FeaturePlot(endo_fibo_object, feature=schwann_markers$gene[1:4])
FeaturePlot(endo_fibo_object, feature=schwann_markers$gene[5:8])
DotPlot(endo_fibo_object, features = schwann_markers$gene[1:20])

NRXN1_high<-endo_fibo_object@meta.data[which(endo_fibo_object$seurat_cluster.008um==3),]

object$final_anno[which(colnames(object)%in%rownames(NRXN1_high))]<-"NRXN1+ HSC"




# Prepare spot coordinates
coords <- as.data.frame(cbind(object@meta.data, object@images$slice1$centroids@coords))
coords$bin<-rownames(coords)

##add NRXN1 expression
gene_exp<-FetchData(object, vars="NRXN1")
gene_exp$cell<-rownames(gene_exp)

rm(object)
gc()

coords<-merge(gene_exp, coords, by.x="cell", by.y="bin")
gc()


 
bin_size <- 13.4543
half_bin <- bin_size / 2
unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 1000
scale_length_plot <- scale_length_um * unit_per_um
x_start <- (-max(coords$x))
y_start <- max(coords$y)

get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

## coordiantes
xmin<-6920
xmax<-7900
ymin<-10160
ymax<-12480

          # plot_bins<-ggplot() +
          #   geom_point(data = coords[-grep("NRXN1",coords$final_anno),],aes(x = y, y = -x, color = final_anno),size = 0.1) +
          #   geom_rect(aes(xmin = min(coords$y), xmax = max(coords$y), ymin = (-min(coords$x)), ymax = (-max(coords$x))), fill = "white", alpha=0.5,linewidth = 1) +  
          #   geom_point(data = coords[grep("NRXN1",coords$final_anno),],aes(x = y, y = -x, color = final_anno),size = 0.1) +
          #   annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start+ 300,yend = x_start+300,
          #            color = "black",linewidth = 1) +
          #   annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start + 200,label = paste0(scale_length_um, " µm"),
          #            color = "black",size = 3,hjust = 0.5) +
          #   # Black border around image
          #   geom_rect(aes(xmin = min(coords$y), xmax = max(coords$y), ymin = (-min(coords$x)), ymax = (-max(coords$x))),
          #             fill = NA,color = "black",linewidth = 1) +
          #   coord_fixed() +colscale_cellType+
          #   # zoom area
          #   geom_rect(aes(xmin = ymin, xmax = ymax, ymin = -xmin, ymax = -xmax),
          #             fill = NA,
          #             color = "black",
          #             linewidth = 0.75) +
          #   theme_void()
          # 
          # full_tissue<-plot_grid(plot_bins+theme(legend.position = "none"), get_leg(plot_bins+guides(color = guide_legend(ncol=4,override.aes = list(size = 2)))), ncol=1, rel_heights = c(3,1))
          # full_tissue




plot_bins<-ggplot() +
  geom_point(data = coords,aes(x = y, y = -x), color = "black",size = 0.8) +
  geom_point(data = coords,aes(x = y, y = -x), color = "grey90",size = 0.1) +
  geom_point(data = coords[grep("hepatocyte",coords$final_anno),],aes(x = y, y = -x, color = final_anno),size = 0.1, alpha=0.1) +
  geom_point(data = coords[grep("NRXN1",coords$final_anno),],aes(x = y, y = -x, color = final_anno),size = 0.1) +
  geom_rect(aes(xmin = ymin, xmax = ymax, ymin = -xmin, ymax = -xmax),fill = NA,color = "black",linewidth = 0.75) +  
  annotate("segment",x = y_start-scale_length_plot-100,xend = y_start-100 ,y = x_start+ 300,yend = x_start+300,color = "black",linewidth = 1) +
  annotate("text",x = y_start - (scale_length_plot / 2)-100,y = x_start + 200,label = paste0(scale_length_um, " µm"),color = "black",size = 3,hjust = 0.5) +
  coord_fixed() +colscale_cellType+
  theme_void()
plot_bins

full_tissue<-plot_grid(plot_bins+theme(legend.position = "none"), get_leg(plot_bins+guides(color = guide_legend(ncol=4,override.aes = list(size = 2)))), ncol=1, rel_heights = c(3,1))
full_tissue

save_plts(full_tissue, "HLiCA_C107_fulltissue_tisse_NRXN1", w=8, h=8)




#######################
## zoom with outlines
#######################

## coordiantes
xmin<-6920
xmax<-7900
ymin<-10160
ymax<-12480

### bin annotation data
coords<-coords[which(coords$x>xmin & coords$x<xmax & coords$y>ymin & coords$y<ymax),]
gc()


# load highest res image
tiff_path<-here("/media/redgar/Seagate Portable Drive/visiumHD_liver/high_res_HnE/H1-NC8X8J2_D1.tif")
tiff_res <- readTIFF(tiff_path, as.is = TRUE)
dim(tiff_res)

## from seurat import above
x_range <- xmin:xmax
y_range <- ymin:ymax



# Step 3: Subset the image
# TIFFs in R are typically in [row, col, channel] format
# Note: y is rows (height), x is columns (width)
zoomed_img <- tiff_res[x_range, y_range, , drop = FALSE]  # if it's RGB
dim(zoomed_img)

# Create RGB array to raster
rgb_img <- rgb(zoomed_img[,,1],
               zoomed_img[,,2],
               zoomed_img[,,3])

# Convert to matrix with same dimensions
dim(rgb_img) <- dim(zoomed_img)[1:2]

# Step 5: Plot
plot(1, type = "n", xlim = c(0, dim(rgb_img)[2]), ylim = c(0, dim(rgb_img)[1]), xlab = "", ylab = "", axes = FALSE)
rasterImage(rgb_img, 0, 0, dim(rgb_img)[2], dim(rgb_img)[1])

rm(tiff_res)
gc()


# Get image dimensions
h <- dim(rgb_img)[1]
w <- dim(rgb_img)[2]

# Convert matrix to data frame
df <- melt(rgb_img)
colnames(df) <- c("y", "x", "color")
df$y <- h - df$y + 1  # Flip y-axis


# Flip y for coords too
coords$y_centred <- (coords$y - ymin)+1
coords$x_centred <- (xmin - (coords$x)+1)+(xmax-xmin)


bin_size <- 13.4543
half_bin <- bin_size / 2

coords_rect <- coords |> 
  mutate(
    xmin = y_centred - half_bin,
    xmax = y_centred + half_bin,
    ymin = x_centred - half_bin,
    ymax = x_centred + half_bin
  )

unit_per_um <- bin_size / 8

# Choose scale bar size in μm
scale_length_um <- 100
scale_length_plot <- scale_length_um * unit_per_um

x_start <- 10
y_start <- -10

image_bounds <- data.frame(
  xmin = min(df$x),
  xmax = max(df$x),
  ymin = min(df$y),
  ymax = max(df$y)
)








aling_views<-plot_grid(
  ggplot() +
    # Raster image
    geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
    scale_fill_identity() +
    # # Transparent bin rectangles
    # geom_rect(data = coords_rect,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #           color = alpha("white", 0.3),fill = NA,linewidth = 0.1) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,
             color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),
             color = "black", size = 2, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA,color = "black",linewidth = 1) +
    # PV and CV
    annotate("text", x = 470, y = 550, label = "CV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    annotate("text", x = 2020, y = 490, label = "PV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 460, y0 = 540, a = 280, b = 120, angle = -45 ), color="black") +
    geom_ellipse(aes(x0 = 2025, y0 = 500, a = 70, b = 70, angle = 0 ), color="black") +
    coord_fixed() +
    colscale_cellType +
    theme_void(),
  
  ggplot() +
    geom_point(data = coords,aes(x = y_centred, y = x_centred, color = final_anno),size = 1, shape=15) +
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", alpha=0.5,color = "black",linewidth = 1) +
    geom_point(data = coords[which(coords$final_anno=="NRXN1+ HSC"),],aes(x = y_centred, y = x_centred),color="black",size = 1, shape=15) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10, color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"), color = "black", size = 2, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA,color = "black",linewidth = 1) +
    # PV and CV
    annotate("text", x = 470, y = 550, label = "CV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    annotate("text", x = 2020, y = 490, label = "PV", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 460, y0 = 540, a = 280, b = 120, angle = -45 ), color="black") +
    geom_ellipse(aes(x0 = 2025, y0 = 500, a = 70, b = 70, angle = 0 ), color="black") +
    coord_fixed() +  colscale_cellType +  theme_void()+theme(legend.position = "none"),
  ncol=1, align="v")

aling_views

save_plts(aling_views, "HLiCA_C107_PV_CV_tisse_NRXN1", w=8, h=8)





aling_views_hep_only<-plot_grid(
  ggplot() +
    geom_point(data = coords,aes(x = y_centred, y = x_centred),size = 1, shape=15, color="grey90") +
    geom_point(data = coords[grep("hepatocyte",coords$final_anno),],aes(x = y_centred, y = x_centred, color = final_anno),size = 1, shape=15) +
    #geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", alpha=0.5,color = "black",linewidth = 1) +
    geom_point(data = coords[which(coords$final_anno=="NRXN1+ HSC"),],aes(x = y_centred, y = x_centred),color="black",size = 1, shape=15) +
    # # Scale bar
    # annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10, color = "black", linewidth = 0.5) +
    # annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"), color = "black", size = 2, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA,color = "black",linewidth = 1) +
    # PV and CV
    annotate("text", x = 470, y = 550, label = "PC", fontface=2, color = "black", size = 5, hjust = 0.5) +
    annotate("text", x = 2020, y = 490, label = "PP", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 460, y0 = 540, a = 280, b = 120, angle = -45 ), color="black") +
    geom_ellipse(aes(x0 = 2025, y0 = 500, a = 70, b = 70, angle = 0 ), color="black") +
    coord_fixed() +  colscale_cellType +  theme_void()+theme(legend.position = "none"),
  
  ggplot() +
    # Raster image
    geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE, alpha=1) +
    scale_fill_identity() +
    # # Transparent bin rectangles
    # geom_rect(data = coords_rect,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #           color = alpha("white", 0.3),fill = NA,linewidth = 0.1) +
    # Scale bar
    annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10,
             color = "black", linewidth = 0.5) +
    annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"),
             color = "black", size = 2, hjust = 0.5) +
    # Black border around image
    geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA,color = "black",linewidth = 1) +
    # PV and CV
    annotate("text", x = 470, y = 550, label = "PC", fontface=2, color = "black", size = 5, hjust = 0.5) +
    annotate("text", x = 2020, y = 490, label = "PP", fontface=2, color = "black", size = 5, hjust = 0.5) +
    geom_ellipse(aes(x0 = 460, y0 = 540, a = 280, b = 120, angle = -45 ), color="black") +
    geom_ellipse(aes(x0 = 2025, y0 = 500, a = 70, b = 70, angle = 0 ), color="black") +
    coord_fixed() +
    colscale_cellType +
    theme_void(),
    ncol=1, align="v")

aling_views_hep_only

save_plts(aling_views_hep_only, "HLiCA_C107_PV_CV_tisse_NRXN1_and hepatocytes", w=8, h=8)




ggplot() +
  geom_point(data = coords,aes(x = y_centred, y = x_centred),size = 1, shape=15, color="grey90") +
  geom_point(data = coords[grep("hepatocyte",coords$final_anno),],aes(x = y_centred, y = x_centred, color = final_anno),size = 1, shape=15) +
  #geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", alpha=0.5,color = "black",linewidth = 1) +
  geom_point(data = coords[which(coords$final_anno=="NRXN1+ HSC"),],aes(x = y_centred, y = x_centred),color="black",size = 1, shape=15) +
  # Scale bar
  annotate("segment", x = x_start, xend = x_start + scale_length_plot, y = y_start-10, yend = y_start-10, color = "black", linewidth = 0.5) +
  annotate("text", x = x_start + scale_length_plot / 2, y = y_start - 30, label = paste0(scale_length_um, " µm"), color = "black", size = 2, hjust = 0.5) +
  # Black border around image
  geom_rect(data = image_bounds,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA,color = "black",linewidth = 1) +
  # PV and CV
  annotate("text", x = 470, y = 550, label = "CV", fontface=2, color = "black", size = 5, hjust = 0.5) +
  annotate("text", x = 2020, y = 490, label = "PV", fontface=2, color = "black", size = 5, hjust = 0.5) +
  geom_ellipse(aes(x0 = 460, y0 = 540, a = 280, b = 120, angle = -45 ), color="black") +
  geom_ellipse(aes(x0 = 2025, y0 = 500, a = 70, b = 70, angle = 0 ), color="black") +
  coord_fixed() +  colscale_cellType +  theme_void()+theme(legend.position = "none")
