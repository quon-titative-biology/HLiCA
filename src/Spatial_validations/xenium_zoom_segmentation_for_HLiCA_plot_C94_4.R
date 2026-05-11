## Load Libraries
library(here)
library(Seurat)

library(SCINA)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

library(sp)
library(scales)
library(viridis)
library(ggridges)

library(tiff)

# library(raster)
# # 
# # Sys.setenv(LIBARROW_MINIMAL = "false") 
# # Sys.setenv(ARROW_WITH_ZSTD = "ON")
# # install.packages("arrow", repos = "https://arrow-r-nightly.s3.amazonaws.com")
# library(arrow)
#library(rgeos)


source("../xenium_liver/scripts/00_pretty_plots.R")
source("../xenium_liver/scripts/00_long_functions.R")
source(here("../xenium_liver/scripts/00_pretty_plots.R"))

mute_colors <- function(cols, amount = 0.4) {
  rgb <- grDevices::col2rgb(cols)
  grDevices::rgb(
    rgb[1, ] * (1 - amount) + 255 * amount,
    rgb[2, ] * (1 - amount) + 255 * amount,
    rgb[3, ] * (1 - amount) + 255 * amount,
    maxColorValue = 255
  )
}

myColors_celltype <- c("#7598bd","#b6cdd9","#aed1ce",  
                       "#FFD700","#f0efa8","#f0efa8",
                       "#a63603","#f03b20","#fcbba1",
                       "#0a15f2","#3e44b8","#5d64f5","#929bda","#5d64f5",
                       "#058205","#8e11bf","#67038f",           
                       "#d1100d","#7a4870", "#27751e",  
                       "#8ede85", "#509418",
                       "#ad235c","#f781b0","#ed0ce2","#ad235c","grey",
                       "hotpink","green","#e0a8ce",
                       "#67038f")  

myColors_celltype_muted<-mute_colors(myColors_celltype)

color_possibilities_celltype<-c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Hepatocyte (Cycling)",  
                                "Cholangiocytes","Cholangiocyte (Biliary)","Cholangiocytes (Biliary tree)",
                                "LSEC","LSEC II","VEC",
                                "HSC (Periportal)","HSC (Quiescent)","HSC","HSC (Activated)","Mesenchyme Low Gene Count (Biliary tree)",
                                "NK-like cells","Mature B-cells","Plasma Cells",           
                                "Erythrocytes","Neutrophil", "gd T-cells",  
                                "Cycling (T Cell)", "CD3+ T-cells",
                                "Macrophage MHCII High","KC Like","Mono-Mac","Myeloid","Doublet",
                                "Check Posistion Hepatocyte (Periportal)","Check Posistion Hepatocyte (middle)","Cycling Myeloid",
                                "NRXN1+ stromal cell")  
names(myColors_celltype) <- color_possibilities_celltype
fillscale_cellType <- scale_fill_manual(name="Cell Type",values = myColors_celltype, drop = T, limits=force)
colscale_cellType <- scale_color_manual(name="Cell Type",values = myColors_celltype, drop = T, limits=force)

names(myColors_celltype_muted) <- color_possibilities_celltype
fillscale_cellType_muted <- scale_fill_manual(name="Cell Type",values = myColors_celltype_muted, drop = T, limits=force)
colscale_cellType_muted <- scale_color_manual(name="Cell Type",values = myColors_celltype_muted, drop = T, limits=force)



############################################################################################################################################################################################
load(here("../HLiCA/data","HLiCA_NRXN1_HLiCA.RData"))

NRXN1_meta$Cell<-sapply(1:nrow(NRXN1_meta), function(x){
  paste(NRXN1_meta$sample[x], strsplit(rownames(NRXN1_meta)[x],"_")[[1]][1], sep="_")
})

load(here("../xenium_liver/data","zonation_allcells.RData"))
zonation_scores<-do.call(rbind, zonation_scores)

zonation_scores$Cell<-sapply(1:nrow(zonation_scores), function(x){
  paste(zonation_scores$sample[x], zonation_scores$cell[x], sep="_")
})

zonation_scores$CellType[which(zonation_scores$Cell%in%NRXN1_meta$Cell)]<-"NRXN1+ stromal cell"
zonation_scores_adult<-zonation_scores[which(zonation_scores$sample%in%c("C94_3","C94_4","C95","C101")),]

df<-zonation_scores_adult[which(zonation_scores_adult$sample=="C94_4"),]




###############
### zoomed plot
###############
### add scale bar
xenium_scale_bar_flip<-function(dataframe, size=4){
  list(annotate("rect", xmin=-max(dataframe$centroid_x)-500, xmax=-max(dataframe$centroid_x),
                ymin=min(-dataframe$centroid_y), ymax=(min(-dataframe$centroid_y))*0.995),
       annotate("text", x=-max(dataframe$centroid_x)-250, y=(min(-dataframe$centroid_y))*1.01, label="0.5mm", size=size))
}


xmin<-800
xmax<-1800
ymin<-2000
ymax<-3000

C94_4_NRXN1_centroids<-ggplot() +
  geom_point(aes(-centroid_x,-centroid_y),df, color="black", size=1.5, shape=19)+
  geom_point(aes(-centroid_x,-centroid_y),df, color="grey95", size=0.6, shape=19)+
  geom_point(aes(-centroid_x,-centroid_y, color=CellType),df[which(df$CellType%in%c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Cholangiocytes")),], size=0.2, shape=19)+
  geom_point(aes(-centroid_x, -centroid_y, color = CellType), df[which(df$CellType%in%c("NRXN1+ stromal cell")),], size=0.8, shape=19)+
  colscale_cellType+xenium_scale_bar_flip(df, 2)+theme(legend.position = "none")+
  geom_rect( aes(xmin = -xmin, xmax = -xmax, ymin = -ymin, ymax = -ymax),    fill = NA, color = "black", linewidth = 0.7)  +theme_void() 
C94_4_NRXN1_centroids
save_plts(C94_4_NRXN1_centroids, "C94_4_NRXN1_centroids", w=8, h=4)


df_zoom <- df[which(df$centroid_x>xmin & df$centroid_x<xmax & df$centroid_y>ymin & df$centroid_y<ymax),]

ggplot() +
  geom_point(aes(-centroid_x,-centroid_y),df_zoom, color="black", size=1.5, shape=19)+
  geom_point(aes(-centroid_x,-centroid_y),df_zoom, color="grey95", size=1, shape=19)+
  geom_point(aes(-centroid_x,-centroid_y, color=CellType),df_zoom[which(df_zoom$CellType%in%c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Cholangiocytes")),], size=1, shape=19)+
  geom_point(aes(-centroid_x, -centroid_y, color = CellType), df_zoom[which(df_zoom$CellType%in%c("NRXN1+ stromal cell")),], size=1, shape=19)+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  #theme_void() + 
  colscale_cellType+xenium_scale_bar_flip(df_zoom, 2)+theme(legend.position = "none")


## other 3 for supp
df_other<-zonation_scores_adult[which(zonation_scores_adult$sample!="C94_4"),]

width_df<-df_other %>%
  group_by(sample) %>%
  summarise(centroid_x_range = max(centroid_x, na.rm = TRUE) - min(centroid_x, na.rm = TRUE)  )

plot_NRXN1_list <- lapply(unique(df_other$sample), function(smp) {
  df <- df_other %>% filter(sample == smp)
  ggplot() +
    geom_point(aes(-centroid_x,-centroid_y),df, color="black", size=1.5, shape=19)+
    geom_point(aes(-centroid_x,-centroid_y),df, color="grey95", size=0.6, shape=19)+
    geom_point(aes(-centroid_x,-centroid_y, color=CellType),df[which(df$CellType%in%c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Cholangiocytes")),], size=0.2, shape=19)+
    geom_point(aes(-centroid_x, -centroid_y, color = CellType), df[which(df$CellType%in%c("NRXN1+ stromal cell")),], size=0.8, shape=19)+
    colscale_cellType+xenium_scale_bar_flip(df, 2)+theme_void()+theme(legend.position = "none")
})

combined <- plot_grid(plotlist = plot_NRXN1_list, ncol = length(plot_NRXN1_list), rel_widths = width_df$centroid_x_range)
combined
save_plts(combined, "Other_xenium_NRXN1_centroids", w=16, h=4)


###################
## C94_4
###################
load(here("../xenium_liver/data/C94_4_segmentation_for_plotting_zoom.RData"))
### Plots
plt_shape$color<-as.factor(plt_shape$cell_id)
levels(plt_shape$color)<-sample(brewer.pal(name="Set1", n=9), length(levels(plt_shape$color)), replace=T)


load(here("data/C94_4_centroid_cellSPA_metrics.RData"))
metrics_all_samples$CellType[which(metrics_all_samples$CellType=="LSEC (Periportal)")]<-"VEC"


plt_shape$Cell_order<-seq(1:nrow(plt_shape))

plt_shape_labels<-merge(plt_shape, metrics_all_samples, by.x="cell_id", by.y="cell")
plt_shape_labels<-plt_shape_labels[order(plt_shape_labels$Cell_order),]

plt_shape_labels$neural<-NA


#### NRXN1 cells
load(here("../HLiCA/data","HLiCA_NRXN1_HLiCA.RData"))

C94_4_NRXN1_meta<-NRXN1_meta[which(NRXN1_meta$sample=="C94_4"),]

C94_4_NRXN1_meta$Cell<-sapply(1:nrow(C94_4_NRXN1_meta), function(x){
  strsplit(rownames(C94_4_NRXN1_meta)[x],"_")[[1]][1]
})


plt_shape_labels$CellType[which(plt_shape_labels$cell_id%in%C94_4_NRXN1_meta$Cell)]<-"NRXN1+ stromal cell"

xmin<-800
xmax<-1800
ymin<-2000
ymax<-3000

zoomed_plt<-plt_shape_labels[which(plt_shape_labels$coord_y > ymin  & plt_shape_labels$coord_y < ymax & plt_shape_labels$coord_x > xmin & plt_shape_labels$coord_x < xmax),]



segementation_plot<-ggplot()+
  geom_polygon(aes(-coord_x, -coord_y, group=cell_id),zoomed_plt, fill="grey90",color="grey90")+
  geom_polygon(aes(-coord_x, -coord_y,  group=cell_id,fill=CellType, color=CellType),
               zoomed_plt[which(zoomed_plt$CellType%in%c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Cholangiocytes")),], 
               linewidth=0.5)+
  geom_polygon(aes(-coord_x, -coord_y,  group=cell_id,fill=CellType),
               zoomed_plt[which(zoomed_plt$CellType%in%c("NRXN1+ stromal cell")),], 
               color="black", linewidth=0.5)+
  annotate("rect",xmin = -xmin,xmax = -xmax ,ymin = -ymin,ymax = -ymax,color = "black", fill=NA, linewidth=1) +
  colscale_cellType+fillscale_cellType+theme_void()+
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_rect(fill = NA, color = "white"))
segementation_plot
save_plts(segementation_plot, "C94_4_NRXN1_zoom", w=6, h=4)

segementation_plot_align<-segementation_plot+theme(legend.position = "none")
segementation_plot_align


##in legend call grey other cells


##############
## dapi
##############
library(raster)
bidcell_dapi <- raster("/media/redgar/Seagate\ Portable\ Drive/liver_BIDCell_output/C94_4/dapi_resized.tif")

bidcell_dapi_df_all <-as.data.frame(bidcell_dapi, xy = TRUE) %>% na.omit()
max(bidcell_dapi_df_all$y)

# Dapi adjust
ymin_dapi<-max(bidcell_dapi_df_all$y)-ymin
ymax_dapi<-max(bidcell_dapi_df_all$y)-ymax

cropbox2 <-c(xmin,xmax,ymax_dapi,ymin_dapi) #(xmin, xmax, ymin, ymax)

DEMcrop2 <- crop(bidcell_dapi, cropbox2)

bidcell_dapi_df <-
  as.data.frame(DEMcrop2, xy = TRUE) %>%
  #--- remove cells with NA for any of the layers ---#
  na.omit(
  )


dapi_plt<-ggplot(data = bidcell_dapi_df) +
  geom_raster(aes(x = x, y = y, fill = dapi_resized)) +
  scale_fill_distiller(type = "seq", direction = 1, palette = "Blues")+
  theme_presentation()+theme(legend.position="none",
                             axis.line = element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_blank(),
                             panel.border = element_rect(fill = NA, color = "white"))
dapi_plt



###############
## together
###############

dapi_threshold<-ggplot()+
  geom_polygon(aes(coord_x, -coord_y, group=cell_id),zoomed_plt, fill="grey90",color="grey90")+
  geom_polygon(aes(coord_x, -coord_y,  group=cell_id, color=CellType),
               zoomed_plt[which(zoomed_plt$CellType%in%c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Cholangiocytes","NRXN1+ stromal cell")),], 
               linewidth=0.5)+
  annotate("rect",xmin = xmin,xmax = xmax ,ymin = -ymin,ymax = -ymax,color = "black", fill=NA) +
  colscale_cellType+theme_void()+
  geom_raster(aes(x = x, y =  (y-max(bidcell_dapi_df_all$y)), fill="red"),data = bidcell_dapi_df[which(bidcell_dapi_df$dapi_resized>1500),]) +
  scale_fill_manual(values = c("blue"))+theme(legend.position = "none")

dapi_threshold_only<-ggplot()+
  annotate("rect",xmin = xmin,xmax = xmax ,ymin = -ymin,ymax = -ymax,color = "black", fill=NA) +
  colscale_cellType+theme_void()+
  geom_raster(aes(x = x, y =  (y-max(bidcell_dapi_df_all$y)), fill="red"),data = bidcell_dapi_df[which(bidcell_dapi_df$dapi_resized>1500),]) +
  scale_fill_manual(values = c("blue"))+theme(legend.position = "none")



aligned<-plot_grid(segementation_plot_align, dapi_plt,dapi_threshold_only,dapi_threshold, ncol=1)
save_plts(aligned, "C94_4_NRXN1_zoom_aligned", w=8, h=20)




##############################
# load H&E
##############################

#tiff_path<-here("/media/redgar/Seagate Portable Drive/xenium_liver/H&E/C94-4_RAW_ch00.tif")
tiff_path<-here("../../Downloads/C94-4_RAW_ch00_rot.tif")

tiff_res <- readTIFF(tiff_path, as.is = TRUE)
dim(tiff_res)

#pixel-to-unit ratio from imageJ: 0.3603520

## coordiantes
xmin<-800/0.3603520
xmax<-1800/0.3603520
ymin<-2000/0.3603520
ymax<-3000/0.3603520

## from seurat import above
x_range <- xmin:xmax
y_range <- ymin:ymax

# Step 3: Subset the image
# TIFFs in R are typically in [row, col, channel] format
# Note: y is rows (height), x is columns (width)
zoomed_img <- tiff_res[y_range, x_range, , drop = FALSE]  # if it's RGB
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



# # Get image dimensions
h <- dim(rgb_img)[1]
w <- dim(rgb_img)[2]

# Convert matrix to data frame
df <- melt(rgb_img)
colnames(df) <- c("y", "x", "color")

# # # H&E adjust
# # ymin_dapi<-max(bidcell_dapi_df_all$y)-ymin
# # ymax_dapi<-max(bidcell_dapi_df_all$y)-ymax
# 
# 
# df$adjust_x<-df$x*(998/2776)


# Rescale to range 1 - 2776
x <- zoomed_plt$coord_x
zoomed_plt$x_scaled <- 1 + (x - min(x)) / (max(x) - min(x)) * (h - 1)
y <- zoomed_plt$coord_y
zoomed_plt$y_scaled <- 1 + (y - min(y)) / (max(y) - min(y)) * (h - 1)



ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE) +
  scale_fill_identity() +
  theme(legend.position = "none")+
  geom_polygon(aes((h-x_scaled)*1.3603520, (h-y_scaled)*1.3603520,  group=cell_id),
               zoomed_plt[which(zoomed_plt$CellType%in%c("NRXN1+ stromal cell")),], 
               color="black", linewidth=0.5)

ggplot() +
  theme(legend.position = "none")+
  geom_polygon(aes(h-x_scaled, h-y_scaled,  group=cell_id),
               zoomed_plt[which(zoomed_plt$CellType%in%c("NRXN1+ stromal cell")),], 
               color="black", linewidth=0.5)

range(zoomed_plt$coord_y)
range(df$y)


x <- df_zoom$centroid_x
df_zoom$x_scaled <- 1 + (x - min(x)) / (max(x) - min(x)) * (h - 1)
y <- df_zoom$centroid_y
df_zoom$y_scaled <- 1 + (y - min(y)) / (max(y) - min(y)) * (h - 1)

ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = color), show.legend = FALSE) +
  scale_fill_identity() +
  geom_point(aes(x_scaled,2776-y_scaled),df_zoom, color="black", size=1.5, shape=19)
  


## overlap cell segmentation
ggplot()+
  geom_polygon(aes(-coord_x, -coord_y, group=cell_id),zoomed_plt, fill="grey90",color="grey90")+
  geom_polygon(aes(-coord_x, -coord_y,  group=cell_id,fill=CellType, color=CellType),
               zoomed_plt[which(zoomed_plt$CellType%in%c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Cholangiocytes")),], 
               linewidth=0.5)+
  geom_polygon(aes(-coord_x, -coord_y,  group=cell_id,fill=CellType),
               zoomed_plt[which(zoomed_plt$CellType%in%c("NRXN1+ stromal cell")),], 
               color="black", linewidth=0.5)+
  annotate("rect",xmin = -xmin,xmax = -xmax ,ymin = -ymin,ymax = -ymax,color = "black", fill=NA, linewidth=1) +
  colscale_cellType+fillscale_cellType+theme_void()+
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.border = element_rect(fill = NA, color = "white"))
