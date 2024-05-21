# Load packages and functions for analysis. Henderson Labs package SeuratPipe is used but code elements can be pulled out and isolated for this analysis if needed.
library(Seurat)
source("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/annotation_utils.R")
seu <- readRDS("project/Liver_Cell_Atlas/Data/First_Download/all_harmonized_seurat.rds") # Load data
modules <- readRDS("project/Liver_Cell_Atlas/Analysis_01/00_Atlas_Markers/modules.rds")

dir.create("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results") # Create results dir
# Run module score analysis
seu <- module_score_analysis(seu = seu,modules_group = modules, plot_dir = "project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/", 
                             reduction = "umap", min.cutoff = NA, legend.position = "right", col_pal = NULL, dims_plot = c(1,2), seed = 42,
                             alpha = c(0.1, 0.9), fig.res = 200, pt.size.factor = 1.1, crop = FALSE, raster = FALSE)

dir.create("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/ordered_feature") # Create results dir
# Print feature plots with order= TRUE
module_score_print(seu = seu,modules_group = modules, plot_dir = "project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/ordered_feature/", 
                   reduction = "umap", min.cutoff = NA, legend.position = "right", col_pal = NULL, dims_plot = c(1,2), seed = 42,
                   alpha = c(0.1, 0.9), fig.res = 200, pt.size.factor = 1.1, crop = FALSE, raster = FALSE, order = TRUE)

saveRDS(seu, "project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/seu.rds")





#########################################################################
# Generate marker lists and use CIPR to see if my annotations match that of cluster identity predictor using HPA database

library(Seurat)
library(SeuratPipe)
library(CIPR)
seu <- readRDS("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/seu.rds")

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/cluster_umap.png",width = 9, height = 7, res = 200, units = "in")
dim_plot(seu, group.by = "seurat_clusters", raster = FALSE, label = T)
dev.off()


pdf("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/cluster_highlight.pdf",width = 9, height = 7)
for (i in levels(seu$seurat_clusters)){
  print(DimPlot(seu, raster=FALSE, cells.highlight = colnames(seu)[seu$seurat_clusters == i])+ ggplot2::ggtitle(paste0("Cluster",i)))
}
dev.off()


markers18 <- FindMarkers(seu, logfc.threshold = 0.5, min.pct = 0.25, only.pos = T, ident.1 = 18)

markers <- FindAllMarkers(seu, logfc.threshold = 0.5, min.pct = 0.25, only.pos = T)

write.csv(markers18, "project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/cluster18_markers.csv")
write.csv(markers, "project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/all_cluster_markers.csv")

sig_markers <- markers[markers$p_val_adj < 0.05,]

CIPR(input_dat = sig_markers,
     comp_method = "logfc_dot_product", 
     reference = "hpca", 
     plot_ind = F,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T)

dir.create("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation") # Create results dir
write.csv(CIPR_top_results, "project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/CIPR_top_output.csv")

##################################################################################
# Create Alpha naming csv / metadata
# Current preliminary alpha naming as I need to separate out lineages to fully validate the annotations


library(Seurat)
library(SeuratPipe)
seu <- readRDS("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/seu.rds")

seu$preliminary_annotations <- as.character(seu$seurat_clusters)
seu$preliminary_annotations <- plyr::mapvalues(seu$preliminary_annotations, from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", 
                                               "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36"),
                            to = c("Hepatocyte","T Cell","Mononuclear Phagocyte","Hepatocyte*","Endothelia","Mononuclear Phagocyte","Mesenchyme","Endothelia",
                                   "T Cell", "Innate Lymphoid Cell", "Endothelia","Cholangiocyte/Hepatocyte","Innate Lymphoid Cell","Plasma B Cell","Endothelia",
                                   "B Cell","Mesenchyme","Cycling","Endothelia","Endothelia","Hepatocyte*","Mononuclear Phagocyte","Mesenchyme","Endothelia",
                                   "Hepatocyte* (+mast+RBC)","Hepatocyte (+mast+RBC)","Mononuclear Phagocyte (+Endothelia)", "Plasmacytoid Dendritic Cell (+mast)","Mononuclear Phagocyte",
                                   "Hepatocyte","Hepatocyte","Hepatocyte*","Neutrophil (+RBC)","T Cell","Hepatocyte* (+MP)", "Hepatocyte*","Hepatocyte (+Endothelia)"))




seu$alpha_annotations <- as.character(seu$seurat_clusters)
seu$alpha_annotations <- plyr::mapvalues(seu$alpha_annotations, from = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", 
                                                                                     "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36"),
                                               to = c("Hepatocyte", "Lymphocyte", "Myeloid Cell", "Hepatocyte*", "Endothelia", "Myeloid Cell", "Mesenchyme", "Endothelia",
                                                      "Lymphocyte", "Lymphocyte", "Endothelia", "Epithelia", "Lymphocyte", "Lymphocyte", "Endothelia", "Lymphocyte", "Mesenchyme",
                                                      "Cycling", "Endothelia", "Endothelia", "Hepatocyte*", "Myeloid Cell", "Mesenchyme", "Endothelia", "Hepatocyte*", "Hepatocyte",
                                                      "Myeloid Cell", "Myeloid Cell", "Myeloid Cell", "Hepatocyte", "Hepatocyte", "Hepatocyte*", "Myeloid Cell", "Lymphocyte",
                                                      "Hepatocyte*", "Hepatocyte*", "Hepatocyte"))


df <- data.frame("Preliminary Annotations" =  seu$preliminary_annotations, "Alpha Annotations" = seu$alpha_annotations)
write.csv(df,"project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/Broad_annotations.csv" )
write.csv(df, "project/Liver_Cell_Atlas/Analysis_01/Annotations_1.csv")


#######################################################################################
# Test to see how well the automated methods match my annotations
# I use the raw annotations from each method but also create a new annotation by taking the modal annotation of cells that make a cluster, i.e. the annotation that is
#most dominant for the cells that make a cluster.

library(Seurat)
library(SeuratPipe)
library(gridExtra)
library(grid)
library(gtable)
seu <- readRDS("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/seu.rds")

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}


# Scmap Cell 
seu$scmap_cell_anno_Cluster <- NA
for(i in levels(seu$seurat_clusters)){
  txt <- Modes(seu$scmap_cell_anno[seu$seurat_clusters==i])
  seu$scmap_cell_anno_Cluster[seu$seurat_clusters==i] <- txt
}

table1 <- tableGrob(table(seu$preliminary_annotations, seu$scmap_cell_anno))
table2 <- tableGrob(table(seu$preliminary_annotations, seu$scmap_cell_anno_Cluster))


title <- textGrob("scMap Cell", gp = gpar(fontsize=50))
subtitle <- textGrob("Raw Numbers", x=0, hjust=0, gp=gpar(fontsize=25, fontface="italic"))

table1 <- gtable_add_rows(table1, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table1 <- gtable_add_rows(table1, heights = grobHeight(title) - unit(60,"mm"), pos = 0)
table1 <- gtable_add_grob(table1, list(title, subtitle), t=c(1,2), l=c(1,1), r=ncol(table1))

subtitle <- textGrob("Mode of Annotation", x=0, hjust=0,  gp=gpar(fontsize=25, fontface="italic"))
table2 <- gtable_add_rows(table2, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table2 <- gtable_add_grob(table2, list(subtitle), t=1, l=1, r=ncol(table2))

pdf("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/scmap_cell.pdf", width = 35, height = 20)
grid.arrange(table1, table2)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/scmap_cell_cluster_umap.png", width = 10, height = 7, units = "in", res = 200)
dim_plot(seu, group.by = "scmap_cell_anno_Cluster", raster = FALSE)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/scmap_cell_cluster_umap_split.png", width = 25, height = 30, units = "in", res = 200)
dim_plot(seu, group.by = "scmap_cell_anno_Cluster", split.by = "scmap_cell_anno_Cluster", raster = FALSE, ncol = 3)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/scmap_cell_umap_split.png", width = 25, height = 30, units = "in", res = 200)
dim_plot(seu, group.by = "scmap_cell_anno", split.by = "scmap_cell_anno", raster = FALSE, ncol = 3)
dev.off()




# Scmap Cluster 
seu$scmap_cluster_anno_Cluster <- NA
for(i in levels(seu$seurat_clusters)){
  txt <- Modes(seu$scmap_cluster_anno[seu$seurat_clusters==i])
  seu$scmap_cluster_anno_Cluster[seu$seurat_clusters==i] <- txt
}

table1 <- tableGrob(table(seu$preliminary_annotations, seu$scmap_cluster_anno))
table2 <- tableGrob(table(seu$preliminary_annotations, seu$scmap_cluster_anno_Cluster))


title <- textGrob("scMap Cluster", gp = gpar(fontsize=50))
subtitle <- textGrob("Raw Numbers", x=0, hjust=0, gp=gpar(fontsize=25, fontface="italic"))

table1 <- gtable_add_rows(table1, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table1 <- gtable_add_rows(table1, heights = grobHeight(title) - unit(60,"mm"), pos = 0)
table1 <- gtable_add_grob(table1, list(title, subtitle), t=c(1,2), l=c(1,1), r=ncol(table1))

subtitle <- textGrob("Mode of Annotation", x=0, hjust=0,  gp=gpar(fontsize=25, fontface="italic"))
table2 <- gtable_add_rows(table2, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table2 <- gtable_add_grob(table2, list(subtitle), t=1, l=1, r=ncol(table2))

pdf("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/scmap_cluster.pdf", width = 35, height = 20)
grid.arrange(table1, table2)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/scmap_cluster_cluster_umap.png", width = 10, height = 7, units = "in", res = 200)
dim_plot(seu, group.by = "scmap_cluster_anno_Cluster", raster = FALSE)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/scmap_cluster_cluster_umap_split.png", width = 25, height = 30, units = "in", res = 200)
dim_plot(seu, group.by = "scmap_cluster_anno_Cluster", split.by = "scmap_cluster_anno_Cluster", raster = FALSE, ncol = 3)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/scmap_cluster_umap_split.png", width = 25, height = 30, units = "in", res = 200)
dim_plot(seu, group.by = "scmap_cluster_anno", split.by = "scmap_cluster_anno", raster = FALSE, ncol = 3)
dev.off()






# Condition Lab
seu$consistent_labs_Cluster <- NA
for(i in levels(seu$seurat_clusters)){
  txt <- Modes(seu$consistent_labs[seu$seurat_clusters==i])
  seu$consistent_labs_Cluster[seu$seurat_clusters==i] <- txt
}

table1 <- tableGrob(table(seu$preliminary_annotations, seu$consistent_labs))
table2 <- tableGrob(table(seu$preliminary_annotations, seu$consistent_labs_Cluster))


title <- textGrob("Consistent Labs", gp = gpar(fontsize=50))
subtitle <- textGrob("Raw Numbers", x=0, hjust=0, gp=gpar(fontsize=25, fontface="italic"))

table1 <- gtable_add_rows(table1, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table1 <- gtable_add_rows(table1, heights = grobHeight(title) - unit(60,"mm"), pos = 0)
table1 <- gtable_add_grob(table1, list(title, subtitle), t=c(1,2), l=c(1,1), r=ncol(table1))

subtitle <- textGrob("Mode of Annotation", x=0, hjust=0,  gp=gpar(fontsize=25, fontface="italic"))
table2 <- gtable_add_rows(table2, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table2 <- gtable_add_grob(table2, list(subtitle), t=1, l=1, r=ncol(table2))

pdf("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/consistent_labs.pdf", width = 35, height = 20)
grid.arrange(table1, table2)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/consistent_labs_cluster_umap.png", width = 10, height = 7, units = "in", res = 200)
dim_plot(seu, group.by = "consistent_labs_Cluster", raster = FALSE)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/consistent_labs_cluster_umap_split.png", width = 25, height = 30, units = "in", res = 200)
dim_plot(seu, group.by = "consistent_labs_Cluster", split.by = "consistent_labs_Cluster", raster = FALSE, ncol = 3)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/consistent_labs_umap_split.png", width = 25, height = 35, units = "in", res = 200)
dim_plot(seu, group.by = "consistent_labs", split.by = "consistent_labs", raster = FALSE, ncol = 3)
dev.off()





# General Lab
seu$general_labs_Cluster <- NA
for(i in levels(seu$seurat_clusters)){
  txt <- Modes(seu$general_labs[seu$seurat_clusters==i])
  seu$general_labs_Cluster[seu$seurat_clusters==i] <- txt
}

table1 <- tableGrob(table(seu$preliminary_annotations, seu$general_labs))
table2 <- tableGrob(table(seu$preliminary_annotations, seu$general_labs_Cluster))


title <- textGrob("General Labs", gp = gpar(fontsize=50))
subtitle <- textGrob("Raw Numbers", x=0, hjust=0, gp=gpar(fontsize=25, fontface="italic"))

table1 <- gtable_add_rows(table1, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table1 <- gtable_add_rows(table1, heights = grobHeight(title) - unit(60,"mm"), pos = 0)
table1 <- gtable_add_grob(table1, list(title, subtitle), t=c(1,2), l=c(1,1), r=ncol(table1))

subtitle <- textGrob("Mode of Annotation", x=0, hjust=0,  gp=gpar(fontsize=25, fontface="italic"))
table2 <- gtable_add_rows(table2, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table2 <- gtable_add_grob(table2, list(subtitle), t=1, l=1, r=ncol(table2))

pdf("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/general_labs.pdf", width = 25, height = 20)
grid.arrange(table1, table2)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/general_labs_cluster_umap.png", width = 10, height = 7, units = "in", res = 200)
dim_plot(seu, group.by = "general_labs_Cluster", raster = FALSE)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/general_labs_cluster_umap_split.png", width = 25, height = 30, units = "in", res = 200)
dim_plot(seu, group.by = "general_labs_Cluster", split.by = "general_labs_Cluster", raster = FALSE, ncol = 3)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/general_labs_umap_split.png", width = 25, height = 30, units = "in", res = 200)
dim_plot(seu, group.by = "general_labs", split.by = "general_labs", raster = FALSE, ncol = 3)
dev.off()







# Marker Lab
seu$marker_labs_Cluster <- NA
for(i in levels(seu$seurat_clusters)){
  txt <- Modes(seu$marker_labs[seu$seurat_clusters==i])
  seu$marker_labs_Cluster[seu$seurat_clusters==i] <- txt
}

table1 <- tableGrob(table(seu$preliminary_annotations, seu$marker_labs))
table2 <- tableGrob(table(seu$preliminary_annotations, seu$marker_labs_Cluster))


title <- textGrob("Marker Labs", gp = gpar(fontsize=50))
subtitle <- textGrob("Raw Numbers", x=0, hjust=0, gp=gpar(fontsize=25, fontface="italic"))

table1 <- gtable_add_rows(table1, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table1 <- gtable_add_rows(table1, heights = grobHeight(title) - unit(60,"mm"), pos = 0)
table1 <- gtable_add_grob(table1, list(title, subtitle), t=c(1,2), l=c(1,1), r=ncol(table1))

subtitle <- textGrob("Mode of Annotation", x=0, hjust=0,  gp=gpar(fontsize=25, fontface="italic"))
table2 <- gtable_add_rows(table2, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table2 <- gtable_add_grob(table2, list(subtitle), t=1, l=1, r=ncol(table2))

pdf("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/marker_labs.pdf", width = 25, height = 20)
grid.arrange(table1, table2)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/marker_labs_cluster_umap.png", width = 10, height = 7, units = "in", res = 200)
dim_plot(seu, group.by = "marker_labs_Cluster", raster = FALSE)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/marker_labs_cluster_umap_split.png", width = 25, height = 30, units = "in", res = 200)
dim_plot(seu, group.by = "marker_labs_Cluster", split.by = "marker_labs_Cluster", raster = FALSE, ncol = 3)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/marker_labs_umap_split.png", width = 25, height = 30, units = "in", res = 200)
dim_plot(seu, group.by = "marker_labs", split.by = "marker_labs", raster = FALSE, ncol = 3)
dev.off()






# Marker General Lab
seu$marker_general_labs_Cluster <- NA
for(i in levels(seu$seurat_clusters)){
  txt <- Modes(seu$marker_general_labs[seu$seurat_clusters==i])
  seu$marker_general_labs_Cluster[seu$seurat_clusters==i] <- txt
}

table1 <- tableGrob(table(seu$preliminary_annotations, seu$marker_general_labs))
table2 <- tableGrob(table(seu$preliminary_annotations, seu$marker_general_labs_Cluster))


title <- textGrob("Marker General Labs", gp = gpar(fontsize=50))
subtitle <- textGrob("Raw Numbers", x=0, hjust=0, gp=gpar(fontsize=25, fontface="italic"))

table1 <- gtable_add_rows(table1, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table1 <- gtable_add_rows(table1, heights = grobHeight(title) - unit(60,"mm"), pos = 0)
table1 <- gtable_add_grob(table1, list(title, subtitle), t=c(1,2), l=c(1,1), r=ncol(table1))

subtitle <- textGrob("Mode of Annotation", x=0, hjust=0,  gp=gpar(fontsize=25, fontface="italic"))
table2 <- gtable_add_rows(table2, heights = grobHeight(subtitle) + unit(58,"mm"), pos = 0)
table2 <- gtable_add_grob(table2, list(subtitle), t=1, l=1, r=ncol(table2))

pdf("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/marker_general_labs.pdf", width = 25, height = 20)
grid.arrange(table1, table2)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/marker_general_labs_cluster_umap.png", width = 10, height = 7, units = "in", res = 200)
dim_plot(seu, group.by = "marker_general_labs_Cluster", raster = FALSE)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/marker_general_labs_cluster_umap_split.png", width = 25, height = 30, units = "in", res = 200)
dim_plot(seu, group.by = "marker_general_labs_Cluster", split.by = "marker_general_labs_Cluster", raster = FALSE, ncol = 3)
dev.off()

png("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/marker_general_labs_umap_split.png", width = 25, height = 30, units = "in", res = 200)
dim_plot(seu, group.by = "marker_general_labs", split.by = "marker_general_labs", raster = FALSE, ncol = 3)
dev.off()













# I do not do this mode approach for sctype output. Only found this annotation after looking at the above and this approach seems to work quite well because it 
# defines cells based on clusters than cells doing something similar by taking the the summation of the scores for each cell for each annoation and picks the highest sumed score.

table1 <- tableGrob(table(seu$preliminary_annotations, seu$customclassif))


title <- textGrob("sc-type", gp = gpar(fontsize=50))

table1 <- gtable_add_rows(table1, heights = grobHeight(title) + unit(60,"mm"), pos = 0)
table1 <- gtable_add_grob(table1, list(title), t=1, l=1, r=ncol(table1))


pdf("project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/results/automation_validation/sctype.pdf", width = 25, height = 10)
grid.arrange(table1)
dev.off()


