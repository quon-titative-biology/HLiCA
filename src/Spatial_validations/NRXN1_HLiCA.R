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

options(future.globals.maxSize = 8000 * 1024^2)

myColors_celltype <- c("#2b8cbe", "#7a0177","#7bccc4",  
                       "#cccc16","#f0efa8","#f0efa8",
                       "#a63603","#f03b20","#fcbba1",
                       "#0a15f2","#3e44b8","#5d64f5","#929bda","#5d64f5",
                       "#058205","#8e11bf","#67038f",           
                       "#d1100d","#7a4870", "#27751e",  
                       "#8ede85", "#509418",
                       "#ad235c","#f781b0","#ed0ce2","#ad235c","grey",
                       "hotpink","green","#e0a8ce",
                       "black")  


color_possibilities_celltype<-c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Hepatocyte (Cycling)",  
                                "Cholangiocytes","Cholangiocyte (Biliary)","Cholangiocytes (Biliary tree)",
                                "LSEC","LSEC II","VEC",
                                "HSC (Periportal)","HSC (Quiescent)","HSC","HSC (Activated)","Mesenchyme Low Gene Count (Biliary tree)",
                                "NK-like cells","Mature B-cells","Plasma Cells",           
                                "Erythrocytes","Neutrophil", "gd T-cells",  
                                "Cycling (T Cell)", "CD3+ T-cells",
                                "Macrophage MHCII High","KC Like","Mono-Mac","Myeloid","Doublet",
                                "Check Posistion Hepatocyte (Periportal)","Check Posistion Hepatocyte (middle)","Cycling Myeloid",
                                "NRXN1+ HSC")  
names(myColors_celltype) <- color_possibilities_celltype
fillscale_cellType <- scale_fill_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)
colscale_cellType <- scale_color_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)

load(file=here("../xenium_liver/data/cell_type_labels_BIDCell.RData"))

print("Load xenium BIDCell Data")

count_files<-c("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C105/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C85/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C95/cell_gene_matrices/2024_07_04_12_35_29/expr_mat.csv",
               "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C101/cell_gene_matrices/2024_06_28_09_23_40/expr_mat.csv")
samples<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")
models<-c("2024_03_21_16_22_10","2024_03_21_16_38_48","2024_03_21_17_03_05","2024_05_17_13_32_03","2024_05_17_13_29_40","2024_07_04_12_35_29","2024_06_28_09_23_40")



d10x.list <- sapply(count_files, function(file_path){
  counts<-read.csv(file_path)
  
  counts$X<-NULL
  rownames(counts) <- counts$cell_id
  counts$cell_id <- NULL
  
  # Transpose the data and convert to sparse matrix.
  mat <- as(t(as.matrix(counts)), "sparseMatrix")
  
  seu <- CreateSeuratObject(counts=mat)
  seu$sample<-strsplit(file_path,"/")[[1]][6]
  seu
})

d10x.list

xenium.obj <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "xenium_liver")
rm(d10x.list)
gc()


metadata_add<-xenium.obj@meta.data

## make seurat object for merging from xenium data
subset.matrix <- xenium.obj@assays$RNA
rm(xenium.obj)
gc()

xenium.obj_RNA <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
identical(colnames(xenium.obj_RNA), rownames(metadata_add))
xenium.obj_RNA<-AddMetaData(xenium.obj_RNA, metadata_add)

xenium.obj_RNA<-JoinLayers(xenium.obj_RNA)



xenium.obj_RNA <- NormalizeData(xenium.obj_RNA)
xenium.obj_RNA <- FindVariableFeatures(xenium.obj_RNA, selection.method = "vst", nfeatures = 2000)
xenium.obj_RNA <- ScaleData(xenium.obj_RNA, verbose = FALSE)
xenium.obj_RNA <- RunPCA(xenium.obj_RNA, npcs = 30, verbose = FALSE)
xenium.obj_RNA <- RunUMAP(xenium.obj_RNA, reduction = "pca", dims = 1:30)
xenium.obj_RNA <- FindNeighbors(xenium.obj_RNA, reduction = "pca", dims = 1:30)
xenium.obj_RNA <- FindClusters(xenium.obj_RNA, resolution = 0.5)

DimPlot(xenium.obj_RNA, label=T)
DimPlot(xenium.obj_RNA)+scale_color_manual(values=c(rep("grey",18),"red", rep("grey",4)))

schwann_markers<-read.csv(here("data/DEG/NRXN1+ Hepatic Stellate Cell.csv"))
schwann_markers<-schwann_markers[which(schwann_markers$avg_log2FC>0),]

DotPlot(object = xenium.obj_RNA, features = schwann_markers$gene[1:600])

table(xenium.obj_RNA$seurat_clusters)


### Neuron
neuron<-c("GPM6B","CDH19","PMP22","SAMHD1",
          "CHL1","NRXN1","PLP1","S100B",
          "CRYAB","S100A10","LGI4","SOX10",
          "IFNGR2","MPZ","FXYD1","GFRA3",
          "SCN7A","PHLDA3")
DotPlot(object = xenium.obj_RNA, features = neuron)


###########################
## highest in neuron markers
###########################
xenium.obj_neuron<-subset(xenium.obj_RNA, subset = seurat_clusters %in% c(14,7))

table(xenium.obj_neuron$seurat_clusters)

xenium.obj_neuron <- NormalizeData(xenium.obj_neuron)
xenium.obj_neuron <- FindVariableFeatures(xenium.obj_neuron, selection.method = "vst", nfeatures = 2000)
xenium.obj_neuron <- ScaleData(xenium.obj_neuron)

xenium.obj_neuron <- RunPCA(xenium.obj_neuron, npcs = 30, features = rownames(xenium.obj_neuron))
xenium.obj_neuron <- RunUMAP(xenium.obj_neuron, dims = 1:30)
xenium.obj_neuron <- FindNeighbors(xenium.obj_neuron, reduction = "pca", dims = 1:30)
xenium.obj_neuron <- FindClusters(xenium.obj_neuron, resolution = 0.2)

BIDCell_UMAP<-DimPlot(xenium.obj_neuron,raster=FALSE,label=T)
BIDCell_UMAP


### Neuron
neuron<-c("GPM6B",
          "CDH19",
          "PMP22",
          "SAMHD1",
          "CHL1",
          "NRXN1",
          "PLP1",
          "S100B",
          "CRYAB",
          "S100A10",
          "LGI4",
          "SOX10",
          "IFNGR2",
          "MPZ",
          "FXYD1",
          "GFRA3",
          "SCN7A",
          "PHLDA3")
FeaturePlot(xenium.obj_neuron, features = neuron)
DotPlot(object = xenium.obj_neuron, features = neuron)



FeaturePlot(xenium.obj_neuron, features=c("HES1", "COL1A1", "COL3A1","PDGFRA","AOX1","ADIPOR1"), raster=F)
DotPlot(object = xenium.obj_neuron, features = c("HES1", "COL1A1", "COL3A1","PDGFRA","AOX1","ADIPOR1"))

DotPlot(object = xenium.obj_neuron, features = c("HES1", "COL1A1", "COL3A1","PDGFRA","AOX1","ADIPOR1"), group.by = "reference_cluster")


custom_markers[grep("CYP",custom_markers$Gene),]
DotPlot(object = xenium.obj_neuron, features = c("CYP1A2", "CYP2A6"))
FeaturePlot(xenium.obj_neuron, features=c("CYP2A6", "CYP1A2", "CYP2A7","CYP2E1"), raster=F) 


####
xenium.obj_neuron<-subset(xenium.obj_neuron, subset = seurat_clusters %in% c(7))
table(xenium.obj_neuron$sample)

neural_meta<-xenium.obj_neuron@meta.data
save(neural_meta, file="data/HLiCA_NRXN1_HLiCA.RData")



load(here("data","HLiCA_NRXN1_HLiCA.RData"))

neural_meta$Cell<-sapply(1:nrow(neural_meta), function(x){
  paste(neural_meta$sample[x], strsplit(rownames(neural_meta)[x],"_")[[1]][1], sep="_")
})

load(here("../xenium_liver/data","zonation_allcells.RData"))
zonation_scores<-do.call(rbind, zonation_scores)

zonation_scores$Cell<-sapply(1:nrow(zonation_scores), function(x){
  paste(zonation_scores$sample[x], zonation_scores$cell[x], sep="_")
})

zonation_scores$CellType[which(zonation_scores$Cell%in%neural_meta$Cell)]<-"NRXN1+ HSC"



ggplot()+
  geom_point(aes(centroid_x,-centroid_y, color=CellType),zonation_scores,size=0.25, shape=19)+
  geom_point(aes(centroid_x,-centroid_y, color=CellType),zonation_scores[which(zonation_scores$CellType=="NRXN1+ HSC"),],size=0.25, shape=19)+
  theme_void()+facet_wrap(~sample, scales = "free")+colscale_cellType



# 
# 
# ###################
# ## C94_4
# ###################
# load(here("data/C94_4_segmentation_for_plotting_zoom.RData"))
# ### Plots
# plt_shape$color<-as.factor(plt_shape$cell_id)
# levels(plt_shape$color)<-sample(brewer.pal(name="Set1", n=9), length(levels(plt_shape$color)), replace=T)
# 
# 
# load(here("data/C94_4_centroid_cellSPA_metrics.RData"))
# metrics_all_samples$CellType[which(metrics_all_samples$CellType=="LSEC (Periportal)")]<-"VEC"
# 
# 
# plt_shape$Cell_order<-seq(1:nrow(plt_shape))
# 
# plt_shape_labels<-merge(plt_shape, metrics_all_samples, by.x="cell_id", by.y="cell")
# plt_shape_labels<-plt_shape_labels[order(plt_shape_labels$Cell_order),]
# 
# plt_shape_labels$neural<-NA
# neural_meta_C94_4<-neural_meta[which(neural_meta$sample=="C94_4"),]
# neural_meta_C94_4$cell<-sapply(1:length(neural_meta_C94_4$cell), function(x) strsplit(neural_meta_C94_4$cell[x], "_")[[1]][1])
# 
# plt_shape_labels$neural[which(plt_shape_labels$cell_id%in%neural_meta_C94_4$cell)]<-"Neural Like Cluster"
# 
# 
# 
# legend_toscale<-get_leg(ggplot(plt_shape_labels, aes(coord_x, -coord_y, color=CellType,  group=cell_id))+
#                           geom_polygon(aes(fill=CellType))+
#                           theme_presentation()+guides(color=guide_legend(ncol=1),fill=guide_legend(ncol=1))+
#                           colscale_cellType+fillscale_cellType+ theme(legend.position="right",
#                                                                       legend.text = element_text(size=10),
#                                                                       legend.title = element_text(size=10),
#                                                                       legend.key.size = unit(0.25, "cm")))
# 
# neural_C94_4<-plot_grid(
#   ggplot()+
#     geom_polygon(aes(coord_x, -coord_y, group=cell_id),plt_shape_labels, fill="grey30", color="grey10")+
#     geom_polygon(aes(coord_x, -coord_y, color=CellType,  group=cell_id,fill=CellType),plt_shape_labels[which(plt_shape_labels$neural%in%c("Neural Like Cluster")),], alpha=0.5)+
#     theme_presentation()+
#     theme(legend.position="none",
#           axis.line = element_blank(),  
#           axis.text = element_blank(),  
#           axis.ticks = element_blank(),  
#           axis.title = element_blank(),  
#           panel.border = element_rect(fill = NA, color = "white"))+colscale_cellType+fillscale_cellType,
#   legend_toscale, ncol=2, rel_widths  = c(5,1))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #############
# ## BIDCell Segmentation
# #############
# 
# 
# xmin<-4000
# xmax<-5000
# ymin<-2600
# ymax<-3100
# 
# 
# zoomed_plt<-plt_shape_labels[which(plt_shape_labels$coord_y > ymin  & plt_shape_labels$coord_y < ymax & plt_shape_labels$coord_x > xmin & plt_shape_labels$coord_x < xmax),]
# 
# legend_toscale<-get_leg(ggplot(zoomed_plt, aes(coord_x, -coord_y, color=CellType,  group=cell_id))+
#                           geom_polygon(aes(fill=CellType))+
#                           theme_presentation()+guides(color=guide_legend(ncol=1),fill=guide_legend(ncol=1))+
#                           colscale_cellType+fillscale_cellType+ theme(legend.position="right",
#                                                                       legend.text = element_text(size=20),
#                                                                       legend.title = element_text(size=20),
#                                                                       legend.key.size = unit(1, "cm")))
# BID_Cell_zoom<-plot_grid(
#   ggplot()+
#     geom_polygon(aes(coord_x, -coord_y, group=cell_id),zoomed_plt, fill="grey20", color="grey5")+
#     geom_polygon(aes(coord_x, -coord_y, color=CellType,  group=cell_id,fill=CellType),zoomed_plt[which(zoomed_plt$neural%in%c("Neural Like Cluster")),], alpha=0.5)+
#     theme_presentation()+
#     theme(legend.position="none",
#           axis.line = element_blank(),  
#           axis.text = element_blank(),  
#           axis.ticks = element_blank(),  
#           axis.title = element_blank(),  
#           panel.border = element_rect(fill = NA, color = "white"))+colscale_cellType+fillscale_cellType,
#   legend_toscale, ncol=2, rel_widths  = c(5,1))
# BID_Cell_zoom
# save_plts_black(BID_Cell_zoom,"BIDCell_C94_4_zoom_neural", w=25, h=15)
# 
# 
# legend.key.size = unit(1, "cm")))
# BID_Cell_zoom<-plot_grid(
#   ggplot()+
#     geom_polygon(aes(coord_x, -coord_y, color=CellType,  group=cell_id,fill=CellType),zoomed_plt, alpha=0.5)+
#     theme_presentation()+
#     theme(legend.position="none",
#           axis.line = element_blank(),  
#           axis.text = element_blank(),  
#           axis.ticks = element_blank(),  
#           axis.title = element_blank(),  
#           panel.border = element_rect(fill = NA, color = "white"))+colscale_cellType+fillscale_cellType,
#   legend_toscale, ncol=2, rel_widths  = c(5,1))
# BID_Cell_zoom
# save_plts_black(BID_Cell_zoom,"BIDCell_C94_4_zoom_neural_allcell", w=25, h=15)
# 
# 
# 
# 
# 
# ###############
# ## Load 10X 
# ###############
# 
# 
# load("/home/redgar/Documents/xenium_liver/data/C94_4_object_raw.RData")
# 
# 
# cropped.coords <- Crop(xenium.obj[["fov"]], x = c(ymin, ymax), y = c(xmin, xmax), coords = "plot")
# xenium.obj[["zoom"]] <- cropped.coords
# # visualize cropped area with cell segmentations & selected molecules
# DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
# 
# ######################
# ### All transcripts in region
# ######################
# 
# gene<-rownames(xenium.obj)
# 
# plt_gene_coor<-do.call(rbind,lapply(1:length(gene), function(x){
#   print(gene[x])
#   if(length(which(names(xenium.obj@images$zoom$molecules)==gene[x]))>0){
#     plt_coor<-as.data.frame(xenium.obj@images$zoom$molecules[which(names(xenium.obj@images$zoom$molecules)==gene[x])][[1]]@coords)
#     plt_coor$gene<-gene[x]
#     plt_coor
#   }
# }))
# 
# 
# 
# transcripts<-ggplot()+
#   geom_point(aes(x,-y),plt_gene_coor,color="white",alpha=0.2, shape=20, stroke=0, size=0.1)+
#   theme_presentation()+theme(legend.position="none",
#                              axis.line = element_blank(),  
#                              axis.text = element_blank(),  
#                              axis.ticks = element_blank(),  
#                              axis.title = element_blank(),  
#                              panel.border = element_rect(fill = NA, color = "white"))
# transcripts
# 
# 
# #### neural markers
# ggplot()+
#   geom_point(aes(x,-y),plt_gene_coor,color="white",alpha=0.2, shape=20, stroke=0, size=0.1)+
#   geom_point(aes(x,-y),plt_gene_coor[which(plt_gene_coor$gene=="LGI4"),],color="#ff0065", size=0.5)+
#   geom_point(aes(x,-y),plt_gene_coor[which(plt_gene_coor$gene=="PMP22"),],color="#00ff1e", size=0.5)+
#   theme_presentation()+colscale_cellType+fillscale_cellType+
#   theme(legend.position="bottom",
#         axis.line = element_blank(),  
#         axis.text = element_blank(),  
#         axis.ticks = element_blank(),  
#         axis.title = element_blank(),  
#         panel.border = element_rect(fill = NA, color = "white"))
# 
# 
# 
# 
# 
# ##############
# ## dapi
# ##############
# bidcell_dapi <- stack("/media/redgar/Seagate\ Portable\ Drive/liver_BIDCell_output/C94_4/dapi_resized.tif")
# 
# bidcell_dapi_df_all <-as.data.frame(bidcell_dapi, xy = TRUE) %>% na.omit()
# max(bidcell_dapi_df_all$y)
# 
# # Dapi adjust
# ymin_dapi<-max(bidcell_dapi_df_all$y)-ymin
# ymax_dapi<-max(bidcell_dapi_df_all$y)-ymax
# 
# cropbox2 <-c(xmin,xmax,ymax_dapi,ymin_dapi) #(xmin, xmax, ymin, ymax)
# 
# DEMcrop2 <- crop(bidcell_dapi, cropbox2)
# 
# bidcell_dapi_df <-
#   as.data.frame(DEMcrop2, xy = TRUE) %>%
#   #--- remove cells with NA for any of the layers ---#
#   na.omit(
#   )
# 
# 
# dapi_plt<-ggplot(data = bidcell_dapi_df) +
#   geom_raster(aes(x = x, y = y, fill = dapi_resized)) +
#   scale_fill_distiller(type = "seq", direction = 1, palette = "Blues")+
#   theme_presentation()+theme(legend.position="none",
#                              axis.line = element_blank(),
#                              axis.text = element_blank(),
#                              axis.ticks = element_blank(),
#                              axis.title = element_blank(),
#                              panel.border = element_rect(fill = NA, color = "white"))
# dapi_plt
# 
# 
# 
# ###############
# ## together
# ###############
# 
# 
# legend_toscale<-get_leg(ggplot()+
#                           geom_point(aes(x,-y),plt_gene_coor,color="white",alpha=0.2, shape=20, stroke=0, size=0.1)+
#                           geom_polygon(aes(coord_x, -coord_y, color=CellType,  group=cell_id),zoomed_plt, fill=NA)+theme_presentation()+
#                           colscale_cellType+fillscale_cellType+
#                           theme_presentation()+theme(legend.position="bottom",
#                                                      legend.text = element_text(size=10),
#                                                      legend.title = element_text(size=10),
#                                                      axis.line = element_blank(),  
#                                                      axis.text = element_blank(),  
#                                                      axis.ticks = element_blank(),  
#                                                      axis.title = element_blank(),  
#                                                      panel.border = element_rect(fill = NA, color = "white")))
# 
# save_plts_black(plot_grid(dapi_plt+ggtitle("Dapi Stain"), legend_toscale,ncol=1, rel_heights = c(5,1)), "C94_4_neural_dapi_zoomed", w=10, h=10)
# save_plts_black(plot_grid(ggplot()+
#                             geom_raster(aes(x = x, y = y, fill="red"),data = bidcell_dapi_df[which(bidcell_dapi_df$dapi_resized>2000),]) +
#                             scale_fill_manual(values = c("blue"))+ggtitle("Dapi Stain Thresholded")+
#                             theme_presentation()+theme(legend.position="none",
#                                                        axis.line = element_blank(),  
#                                                        axis.text = element_blank(),  
#                                                        axis.ticks = element_blank(),  
#                                                        axis.title = element_blank(),  
#                                                        panel.border = element_rect(fill = NA, color = "white")), 
#                           legend_toscale,ncol=1, rel_heights = c(5,1)),
#                 "C94_4_neural_dapi_segment_zoomed", w=10, h=10)
# 
# 
# save_plts_black(  plot_grid(ggplot()+
#                               geom_point(aes(x,-y),plt_gene_coor,color="white",alpha=0.4, shape=20, stroke=0, size=0.5)+
#                               colscale_cellType+fillscale_cellType+ggtitle("All Transcripts")+
#                               theme_presentation()+theme(legend.position="bottom",
#                                                          axis.line = element_blank(),  
#                                                          axis.text = element_blank(),  
#                                                          axis.ticks = element_blank(),  
#                                                          axis.title = element_blank(),  
#                                                          panel.border = element_rect(fill = NA, color = "white")),
#                             legend_toscale,ncol=1, rel_heights = c(5,1)),
#                   "C94_4_neural_zoomed_All_transcripts", w=10, h=10)
# 
# save_plts_black(  plot_grid(ggplot()+
#                               geom_point(aes(x,-y),plt_gene_coor,color="white",alpha=0.4, shape=20, stroke=0, size=0.5)+
#                               geom_point(aes(x,-y),plt_gene_coor[which(plt_gene_coor$gene=="LGI4"),],color="#ff0065", size=1)+
#                               geom_point(aes(x,-y),plt_gene_coor[which(plt_gene_coor$gene=="PMP22"),],color="#00ff1e", size=1)+
#                               colscale_cellType+fillscale_cellType+ggtitle("All Transcripts")+
#                               theme_presentation()+theme(legend.position="bottom",
#                                                          axis.line = element_blank(),  
#                                                          axis.text = element_blank(),  
#                                                          axis.ticks = element_blank(),  
#                                                          axis.title = element_blank(),  
#                                                          panel.border = element_rect(fill = NA, color = "white")),
#                             legend_toscale,ncol=1, rel_heights = c(5,1)),
#                   "C94_4_neural_zoomed_All_transcripts_neural_markers", w=10, h=10)
# 
# 
# 
# 
# 
# 
# save_plts_black(  plot_grid(ggplot()+
#                               geom_point(aes(x,-y),plt_gene_coor,color="white",alpha=0.4, shape=20, stroke=0, size=0.5)+
#                               geom_polygon(aes(coord_x, -coord_y, color=CellType,  group=cell_id),zoomed_plt, linewidth=0.75, fill=NA)+theme_presentation()+
#                               colscale_cellType+fillscale_cellType+ggtitle("All Transcripts BIDCell Segmentation")+
#                               theme_presentation()+theme(legend.position="none",
#                                                          axis.line = element_blank(),  
#                                                          axis.text = element_blank(),  
#                                                          axis.ticks = element_blank(),  
#                                                          axis.title = element_blank(),  
#                                                          panel.border = element_rect(fill = NA, color = "white")),
#                             legend_toscale,ncol=1, rel_heights = c(5,1)),
#                   "C94_4_neural_zoomed_BIDCell", w=10, h=10)
# 
# 
# save_plts_black(plot_grid(ggplot()+
#                             geom_point(aes(x,-y),data=plt_gene_coor, color="white",alpha=0.4, shape=20, stroke=0, size=0.5)+
#                             geom_raster(aes(x = x, y =  (y-max(bidcell_dapi_df_all$y)), fill="red"),data = bidcell_dapi_df[which(bidcell_dapi_df$dapi_resized>1500),]) +
#                             scale_fill_manual(values = c("blue"))+colscale_cellType+ggtitle("DAPI Thresholded BIDCell Segmentation")+
#                             theme_presentation()+theme(legend.position="none",
#                                                        axis.line = element_blank(),  
#                                                        axis.text = element_blank(),  
#                                                        axis.ticks = element_blank(),  
#                                                        axis.title = element_blank(),  
#                                                        panel.border = element_rect(fill = NA, color = "white")),
#                           legend_toscale,ncol=1, rel_heights = c(5,1)),
#                 "C94_4_neural_dapi_trancript_zoomed", w=10, h=10)
# 
# 
# save_plts_black(plot_grid(ggplot()+
#                             geom_point(aes(x,-y),data=plt_gene_coor, color="white",alpha=0.4, shape=20, stroke=0, size=0.5)+
#                             geom_raster(aes(x = x, y =  (y-max(bidcell_dapi_df_all$y)), fill="red"),data = bidcell_dapi_df[which(bidcell_dapi_df$dapi_resized>1500),]) +
#                             geom_polygon(aes(coord_x, -coord_y, color=CellType,  group=cell_id),zoomed_plt, linewidth=0.75, fill=NA)+
#                             scale_fill_manual(values = c("blue"))+colscale_cellType+ggtitle("Dapi plus All Transcripts BIDCell Segmentation")+
#                             theme_presentation()+theme(legend.position="none",
#                                                        axis.line = element_blank(),  
#                                                        axis.text = element_blank(),  
#                                                        axis.ticks = element_blank(),  
#                                                        axis.title = element_blank(),  
#                                                        panel.border = element_rect(fill = NA, color = "white")),
#                           legend_toscale,ncol=1, rel_heights = c(5,1)),
#                 "C94_4_neural_dapi_trancript_zoomed_BIDCell", w=10, h=10)
# 
# 
# 
# save_plts_black(plot_grid(ggplot()+
#                             geom_polygon(aes(coord_x, -coord_y, group=cell_id),zoomed_plt, fill="grey20", color="grey5")+
#                             geom_polygon(aes(coord_x, -coord_y, color=CellType,  group=cell_id),zoomed_plt[which(zoomed_plt$neural=="Neural Like Cluster"),], linewidth=0.75, fill=NA)+
#                             geom_point(aes(x,-y),plt_gene_coor[which(plt_gene_coor$gene=="LGI4"),],color="#ff0065", size=0.5)+
#                             geom_point(aes(x,-y),plt_gene_coor[which(plt_gene_coor$gene=="PMP22"),],color="#00ff1e", size=0.5)+ 
#                             geom_point(aes(x,-y),plt_gene_coor[which(plt_gene_coor$gene=="COL1A1"),],color="cornflowerblue", size=0.5)+
#                             annotate("text",x=4100,y=-3000, label="LGI4", color="#ff0065")+
#                             annotate("text",x=4100,y=-3025, label="PMP22", color="#00ff1e")+
#                             annotate("text",x=4100,y=-3050, label="COL1A1", color="cornflowerblue")+
#                             scale_fill_manual(values = c("blue"))+colscale_cellType+ggtitle("Markers and BIDCell Segmentation")+
#                             theme_presentation()+theme(legend.position="none",
#                                                        axis.line = element_blank(),  
#                                                        axis.text = element_blank(),  
#                                                        axis.ticks = element_blank(),  
#                                                        axis.title = element_blank(),  
#                                                        panel.border = element_rect(fill = NA, color = "white")),
#                           legend_toscale,ncol=1, rel_heights = c(5,1)),
#                 "C94_4_neural_marker_zoomed_segmentation", w=10, h=10)
# 
# 
# 
# save_plts_black(plot_grid(ggplot()+
#                             geom_polygon(aes(coord_x, -coord_y, group=cell_id),zoomed_plt, fill="grey20", color="grey5")+
#                             geom_polygon(aes(coord_x, -coord_y, color=CellType,  group=cell_id,fill=CellType),zoomed_plt[which(zoomed_plt$neural%in%c("Neural Like Cluster")),], alpha=0.5)+
#                             fillscale_cellType+colscale_cellType+ggtitle("BIDCell Segmentation")+
#                             theme_presentation()+theme(legend.position="none",
#                                                        axis.line = element_blank(),  
#                                                        axis.text = element_blank(),  
#                                                        axis.ticks = element_blank(),  
#                                                        axis.title = element_blank(),  
#                                                        panel.border = element_rect(fill = NA, color = "white")),
#                           legend_toscale,ncol=1, rel_heights = c(5,1)),
#                 "C94_4_neural_zoomed_segmentation", w=10, h=10)
# 
# save_plts_black(plot_grid(ggplot()+
#                             geom_polygon(aes(coord_x, -coord_y, color=CellType,  group=cell_id,fill=CellType),zoomed_plt, alpha=0.5)+
#                             fillscale_cellType+colscale_cellType+ggtitle("BIDCell Segmentation")+
#                             theme_presentation()+theme(legend.position="none",
#                                                        axis.line = element_blank(),  
#                                                        axis.text = element_blank(),  
#                                                        axis.ticks = element_blank(),  
#                                                        axis.title = element_blank(),  
#                                                        panel.border = element_rect(fill = NA, color = "white")),
#                           legend_toscale,ncol=1, rel_heights = c(5,1)),
#                 "C94_4_neural_zoomed_segmentation_allcell", w=10, h=10)
