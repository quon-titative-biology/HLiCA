# Load packages and functions for analysis.
library(Seurat)
library(dplyr)
library(tidytext)
library(ggplot2)
library(ggbrace)
library(CIPR)
source("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/annotation_utils.R")
source("project/Functions/Utils.R")

###########################################################################################################
# Load dataset, add in old annotations and cluster the new dataset

seu <- readRDS("project/Liver_Cell_Atlas/Data/Second_Download/healthy_RNA_merged_harmonized.rds") # Load data
analysis_v1_annotations <- read.csv("project/Liver_Cell_Atlas/Analysis_01/Annotations_4.csv", row.names = 1)
seu <- SeuratPipe::add_umap_embedding(seu, "project/Liver_Cell_Atlas/Data/Second_Download/harmony_umap.csv")

# Add in Old annotations and add "New cell" for the new cells added to the dataset and call all the new cells false and potential doublest as the annotation get updated as 
# I go along
colnames(analysis_v1_annotations) <- c("OLD.Preliminary.Annotations", "OLD.Alpha.Annotations", "OLD.Potential.Doublets")
seu <- AddMetaData(seu, analysis_v1_annotations)
seu$OLD.Preliminary.Annotations[is.na(seu$OLD.Preliminary.Annotations)] <- "New Cells"
seu$OLD.Alpha.Annotations[is.na(seu$OLD.Alpha.Annotations)] <- "New Cells"
seu$OLD.Potential.Doublets[is.na(seu$OLD.Potential.Doublets)] <- "FALSE"


table(seu$OLD.Alpha.Annotations)

#Cholangiocyte      Doublets    Endothelia   Erythrocyte    Hepatocyte    Lymphocyte    Mesenchyme  Myeloid Cell     New Cells 
#7954                 31         84401           225        155410        105983         34794         72581        183301 

length(colnames(seu)) - length(rownames(analysis_v1_annotations))
#[1] 183301

# The number of cells with a matching name is the same as the number of cells from the annotations with the same number of cells showing as new cells as the remainder 
# between the number of old annotated cells and new number of cells.

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:40)
seu <- FindClusters(seu, resolution = 0.5) # Using same resolution as was used according to the command log of V1 dataset.
# This gives 37 clusters (similar number as v1 dataset)


dir.create("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results")

plot <- SeuratPipe::dim_plot(seu, group.by = "seurat_clusters", col_pal = SeuratPipe:::discrete_col_pal, raster = FALSE, label=T)
png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Clusters_UMAP.png", width = 8, height = 7, units = "in", res = 300)
print(plot)
dev.off()

saveRDS(seu, "project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/seu.rds")

#################################################################################################################

# Check how well the previous annotations match to this dataset. 
# Any cluster with a 95% match of a single annotation (excluding new cells) will be labelled as that annotation (including new cells) for now.
# This will then all be check at the end and clusters may change annotation

seu <- readRDS("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/seu.rds")

meta.data <- seu@meta.data
meta.data <- meta.data[-which(meta.data$OLD.Alpha.Annotations == "New Cells"),]

col_pal <- c("Myeloid Cell" = "skyblue2","Mesenchyme" = "darkolivegreen4","Lymphocyte"= "tan3","Hepatocyte" = "#FB9A99","Erythrocyte" = "#CAB2D6", "Endothelia" = "#FDBF6F",
             "Doublets" = "#6699CB","Cholangiocyte" = "indianred", "New Cells" = "darkgrey")

# Make empty dataframe to store data used in stacked bar plot
cl_condition <- data.frame(seurat_clusters=factor(),
                           OLD.Alpha.Annotations=character(),
                           n = integer(),
                           freq = numeric(),
                           ordering = integer())
# Fill empty dataframe with data for each cluster. Have to do it cluster by cluster to get the ordering correct so that we have a value to order the labels in each stacked bar
# by percentage.
for (i in levels(meta.data$seurat_clusters)[-which(levels(meta.data$seurat_clusters) == "35")]){
  tmp <- meta.data |> dplyr::group_by(seurat_clusters, OLD.Alpha.Annotations) |>
    dplyr::summarise(n = dplyr::n()) |> dplyr::mutate(freq = n / sum(n)) |> filter(seurat_clusters == i) |> arrange(desc(freq)) |> mutate(ordering = 1:n())
  cl_condition <- rbind(cl_condition, tmp)
}
# Cluster 35 only has new cells and so wont be given a name based on this analysis. Instead I will make a decision of what to call these based on the expression
# of genes from this analysis, or included in any dataset that will have to be seperated out because they have less than 95% of any one particular annotation.

nlevels(droplevels(as.factor(cl_condition$seurat_clusters))) == nlevels(seu$seurat_clusters)-1 # Check I have created all the dfs I need to. Need to remove one for not including cluster 35

# Create a dataframe with all clusters to generate the order in which the clusters should appear with those clusters with closer to 100% of one annotation appearing first.
bar_order_data <- meta.data |> dplyr::group_by(seurat_clusters, OLD.Alpha.Annotations) |>
  dplyr::summarise(n = dplyr::n()) |> dplyr::mutate(freq = n / sum(n)) 

#Change factoring to order clusters so that those with higher percentage annotations appear first
cl_condition$seurat_clusters <- factor(cl_condition$seurat_clusters, levels = unique(as.character(bar_order_data$seurat_clusters[rev(order(bar_order_data$freq))])))

# Generate the ID for ordering within each stacked bar so larger percentage annotations appear first
cl_condition <- cl_condition %>%
  arrange(seurat_clusters, ordering) |>
  mutate(
    id_ordering = paste(seurat_clusters, ordering, sep = "_"),
    id_ordering = forcats::fct_inorder(id_ordering)
  )

plot <- ggplot(cl_condition, aes(
  x = seurat_clusters,
  y = freq,
  fill = OLD.Alpha.Annotations,
  group = id_ordering
)) + theme_classic() +
  geom_col(position = position_stack(reverse = TRUE)) +
  labs(fill = "ordering") +
  guides(fill = guide_legend(reverse = TRUE)) +
  ggplot2::scale_fill_manual(values =col_pal) +
  geom_hline(yintercept = 0.95, color = "black", linewidth=1)



png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Previous_Alpha_Annotation_Frequency_Plot.png", width = 10, height = 7, units = "in", res = 300)
print(plot)
dev.off()

############################################################################################

# repeat of above but including new cells this time

seu <- readRDS("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/seu.rds")

meta.data <- seu@meta.data

col_pal <- c("Myeloid Cell" = "skyblue2","Mesenchyme" = "darkolivegreen4","Lymphocyte"= "tan3","Hepatocyte" = "#FB9A99","Erythrocyte" = "#CAB2D6", "Endothelia" = "#FDBF6F",
             "Doublets" = "#6699CB","Cholangiocyte" = "indianred", "New Cells" = "darkgrey")


# Make empty dataframe to store data used in stacked bar plot
cl_condition <- data.frame(seurat_clusters=factor(),
                           OLD.Alpha.Annotations=character(),
                           n = integer(),
                           freq = numeric(),
                           ordering = integer())
# Fill empty dataframe with data for each cluster. Have to do it cluster by cluster to get the ordering correct so that we have a value to order the labels in each stacked bar
# by percentage.
for (i in levels(meta.data$seurat_clusters)){
  tmp <- meta.data |> dplyr::group_by(seurat_clusters, OLD.Alpha.Annotations) |>
    dplyr::summarise(n = dplyr::n()) |> dplyr::mutate(freq = n / sum(n)) |> filter(seurat_clusters == i) |> arrange(desc(freq)) |> mutate(ordering = 1:n())
  cl_condition <- rbind(cl_condition, tmp)
}

nlevels(droplevels(as.factor(cl_condition$seurat_clusters))) == nlevels(seu$seurat_clusters) # Check I have created all the dfs I need to. 

bar_order_data <- meta.data |> dplyr::group_by(seurat_clusters, OLD.Alpha.Annotations) |>
  dplyr::summarise(n = dplyr::n()) |> dplyr::mutate(freq = n / sum(n))


cl_condition$seurat_clusters <- factor(cl_condition$seurat_clusters, levels = unique(as.character(bar_order_data$seurat_clusters[rev(order(bar_order_data$freq))])))


cl_condition <- cl_condition %>%
  arrange(seurat_clusters, ordering) |>
  mutate(
    id_ordering = paste(seurat_clusters, ordering, sep = "_"),
    id_ordering = forcats::fct_inorder(id_ordering)
  )

plot <- ggplot(cl_condition, aes(
  x = seurat_clusters,
  y = freq,
  fill = OLD.Alpha.Annotations,
  group = id_ordering
)) + theme_classic() +
  geom_col(position = position_stack(reverse = TRUE)) +
  labs(fill = "ordering") +
  guides(fill = guide_legend(reverse = TRUE)) +
  ggplot2::scale_fill_manual(values = col_pal) 

png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Previous_Alpha_Annotation_Frequency_Plot_with_New_Cells.png", width = 10, height = 7, units = "in", res = 300)
print(plot)
dev.off()

####################################################################################################################

# Dot plot of Henderson Gene list. Clustered that have been defined by the annotation threshold defined above are grouped together, with the other clusters put into
# an unknown grouping at this stage. 


seu <- readRDS("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/seu.rds")

col_pal <- c("Myeloid Cell" = "skyblue2","Mesenchyme" = "darkolivegreen4","Lymphocyte"= "tan3","Hepatocyte" = "#FB9A99","Erythrocyte" = "#CAB2D6", "Endothelia" = "#FDBF6F",
             "Doublets" = "#6699CB","Cholangiocyte" = "indianred", "New Cells" = "darkgrey", "Unknown" = "darkgrey")

# Here I list out the temporary naming of each cluster. So any cluster with over 95% single annotation(visualised using the frequency plots generated above) 
# was labelled with that annotation. The rest were labelled as unknown for now 
annotations <- c("Unknown","Hepatocyte","Hepatocyte","Lymphocyte","Myeloid Cell","Unknown","Endothelia","Unknown","Myeloid Cell","Endothelia","Lymphocyte","Endothelia","Unknown",
  "Lymphocyte","Unknown","Lymphocyte","Unknown","Unknown","Unknown","Lymphocyte","Unknown","Endothelia","Unknown","Myeloid Cell","Unknown","Unknown","Hepatocyte",
  "Hepatocyte", "Myeloid Cell","Unknown","Unknown","Lymphocyte","Lymphocyte","Endothelia","Hepatocyte","Unknown","Myeloid Cell")

length(annotations) # check I have the correct number for the amount of clusters
table(annotations) # check all names spelt correctly of each annotation (ie. no single annotation was miss spelt)

# Add these temporaty annotations into the seurat object
seu$Alpha.Annotations.tmp <- plyr::mapvalues(seu$seurat_clusters, from = levels(seu$seurat_clusters), to = annotations)

# Check these annotation match with the frequency plot
tmp <- seu@meta.data
tmp$seurat_clusters <- factor(tmp$seurat_clusters, levels = c("36", "34", "33", "32", "31", "28", "13", "10", "8",  "4",  "3",  "15", "2",  "11", "9",  "6",
                                                              "23", "1",  "21", "19", "26", "27", "5",  "29", "7",  "0",  "12", "25", "24", "14", "17", "30",
                                                              "20", "16", "18", "22","35")) # change order of seurat clusters inline with Annotation frequency plot just to check new annotations easily

table(tmp$seurat_clusters, tmp$Alpha.Annotations.tmp) # check new annotations



# Read In Lineage Gene Signatures
modules <- readRDS("project/Liver_Cell_Atlas/Analysis_02/00_Atlas_Markers/modules.rds")


# Make temp seurat file to create re ordered clusters for dotplot so that those cluster that have adopted a annotation can be grouped by annotation
# With thhose labelled as unknown grouped together at the end. Cluster 35 (made of entirly new cells) was added at the end.
tmp <- seu
tmp$seurat_clusters <- factor(seu$seurat_clusters, levels = c(1,2,26,27,34,3,10,13,15,19,31,32,4,8,23,28,36,6,9,11,21,33,5,29,7,0,12,25,24,14,17,30,20,16,18,22,35))
Idents(tmp) <- "Alpha.Annotations.tmp"

plot <- DotPlot(tmp, features = unique(unlist(modules$Hendo)), group.by = "seurat_clusters") + ggplot2::coord_flip(x=c(1,69), y=c(1,37),clip = "off") + 
  eval(rlang::expr(ggplot2::scale_colour_distiller(name = NULL,type = "div", palette = "RdYlBu"))) + 
  theme(axis.title.x = element_text(margin = margin(t = 80))) + theme(axis.title.y = element_text(margin = margin(r = 120))) +
  ggplot2::theme(legend.position = "right") + 
  geom_brace(xstart = -1, xend = -3, ystart = 1, yend = 5,pointing="side") +
  annotate("text",x=-4,y=3,label="Hepatocytes") +
  geom_brace(xstart = -1, xend = -3, ystart = 6, yend = 12,pointing="side") +
  annotate("text",x=-4,y=9,label="Lymphocyte") +
  geom_brace(xstart = -1, xend = -3, ystart = 13, yend = 17,pointing="side") +
  annotate("text",x=-4,y=15,label="Myeloid Cell") +
  geom_brace(xstart = -1, xend = -3, ystart = 18, yend = 22,pointing="side") +
  annotate("text",x=-4,y=20,label="Endothelia") +
  geom_brace(xstart = -1, xend = -3, ystart = 23, yend = 37,pointing="side") +
  annotate("text",x=-4,y=30,label="Unknown") +

  geom_brace(xstart = 1, xend = 7, ystart = -3, yend = -5) +
  annotate("text",x=3,y=-6.5,label="Hepatocytes",angle = 45) +
  geom_brace(xstart = 8, xend = 14, ystart = -3, yend = -5) +
  annotate("text",x=10.5,y=-6.5,label="Cholangiocytes",angle = 45) +
  geom_brace(xstart = 15, xend = 22, ystart = -3, yend = -5) +
  annotate("text",x=17.5,y=-6.5,label="Endothelia",angle = 45) +
  geom_brace(xstart = 23, xend = 31, ystart = -3, yend = -5) +
  annotate("text",x=26,y=-6.5,label="Mesenchyme",angle = 45) +
  geom_brace(xstart = 32, xend = 33, ystart = -3, yend = -5) +
  annotate("text",x=31.5,y=-6.5,label="Mesothelia",angle = 45) +
  geom_brace(xstart = 34, xend = 37, ystart = -3, yend = -5) +
  annotate("text",x=34.5,y=-6.5,label="Leucocytes",angle = 45) +
  geom_brace(xstart = 38, xend = 43, ystart = -3, yend = -5) +
  annotate("text",x=40.5,y=-5.7,label="MPs",angle = 45) +
  geom_brace(xstart = 44, xend = 46, ystart = -3, yend = -5) +
  annotate("text",x=44,y=-6.5,label="Lymphocytes",angle = 45) +
  geom_brace(xstart = 47, xend = 52, ystart = -3, yend = -5) +
  annotate("text",x=49,y=-6.2,label="T Cells",angle = 45) +
  geom_brace(xstart = 53, xend = 54, ystart = -3, yend = -5) +
  annotate("text",x=53,y=-6.2,label="B Cells",angle = 45) +
  geom_brace(xstart = 55, xend = 56, ystart = -3, yend = -5) +
  annotate("text",x=55.5,y=-5.7,label="ILCs",angle = 45) +
  geom_brace(xstart = 57, xend = 58, ystart = -3, yend = -5) +
  annotate("text",x=57,y=-6.3,label="Mast Cells",angle = 45) +
  geom_brace(xstart = 59, xend = 60, ystart = -3, yend = -5) +
  annotate("text",x=59.5,y=-6,label="PDCs",angle = 45) +
  geom_brace(xstart = 61, xend = 64, ystart = -3, yend = -5) +
  annotate("text",x=61.5,y=-6.5,label="Erythrocytes",angle = 45) +
  geom_brace(xstart = 65, xend = 65, ystart = -3, yend = -5) +
  annotate("text",x=64,y=-6.5,label="Neutrophils",angle = 45) +
  geom_brace(xstart = 66, xend = 69, ystart = -3, yend = -5) +
  annotate("text",x=67,y=-6,label="Cycling",angle = 45) +
  
  theme(plot.margin = unit(c(0.01, 0.01, 0.05, 0.05), units="npc"))

png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Henderson_Module_Gene_by_Cluster_Dotplot.png", width = 15, height = 15, units = "in", res = 300)
print(plot)
dev.off()


plot <- SeuratPipe::subset_dim_plot(seu, subset.by = "Alpha.Annotations.tmp", col_pal = c(col_pal["Unknown"], col_pal["Hepatocyte"],col_pal["Lymphocyte"],
                                                                                         col_pal["Myeloid Cell"],col_pal["Endothelia"]), raster = FALSE) #Have to specify colour ordering as doesnt give correct colours to cell type otherwise
png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Temporary_Alpha_Annotations_UMAP.png", width = 25, height = 15, units = "in", res = 300)
print(plot)
dev.off()


plot <- SeuratPipe::subset_dim_plot(seu, subset.by = "OLD.Alpha.Annotations", col_pal = c(col_pal["Cholangiocyte"], col_pal["Doublets"],col_pal["Endothelia"],
                                                                                      col_pal["Erythrocyte"],col_pal["Hepatocyte"], col_pal["Lymphocyte"],
                                                                                      col_pal["Mesenchyme"],col_pal["Myeloid Cell"], col_pal["New Cells"]), raster = FALSE) #Have to specify colour ordering as doesnt give correct colours to cell type otherwise
png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Previous_Alpha_Annotations_UMAP.png", width = 25, height = 20, units = "in", res = 300)
print(plot)
dev.off()

saveRDS(seu, "project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/seu.rds")

######################################################################################################

# Here I am going to map back all the unknown labelled clusters back to the original dataset to see why clusters in this dataset have a varied annotation spread.
# If they map to many clusters and areas of the previous dataset, that would make sen

seu <- readRDS("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/seu.rds")
old <- readRDS("project/Liver_Cell_Atlas/Data/First_Download/all_harmonized_seurat.rds")

plot1 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "5"], raster = FALSE) + ggtitle("Cluster 5")
plot2 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "29"], raster = FALSE) + ggtitle("Cluster 29")
plot3 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "7"], raster = FALSE) + ggtitle("Cluster 7")
plot4 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "0"], raster = FALSE) + ggtitle("Cluster 0")
plot5 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "12"], raster = FALSE) + ggtitle("Cluster 12")
plot6 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "25"], raster = FALSE) + ggtitle("Cluster 25")
plot7 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "24"], raster = FALSE) + ggtitle("Cluster 24")
plot8 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "14"], raster = FALSE) + ggtitle("Cluster 14")
plot9 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "17"], raster = FALSE) + ggtitle("Cluster 17")
plot10 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "30"], raster = FALSE) + ggtitle("Cluster 30")
plot11 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "20"], raster = FALSE) + ggtitle("Cluster 20")
plot12 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "16"], raster = FALSE) + ggtitle("Cluster 16")
plot13 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "18"], raster = FALSE) + ggtitle("Cluster 18")
plot14 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "22"], raster = FALSE) + ggtitle("Cluster 22")
plot15 <- DimPlot(old, cells.highlight = colnames(seu)[seu$seurat_clusters == "35"], raster = FALSE) + ggtitle("Cluster 35")

png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Map_Unknown_Clusters_to_Previous_Dataset.png", width = 25, height = 20, units = "in", res = 300)
print(patchwork::wrap_plots(list(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, plot13, plot14, plot15), ncol = 4, nrow = 4))
dev.off()


############################################################################################
#  Henderson Labs package SeuratPipe is used but code elements can be pulled out and isolated for this analysis if needed.


seu <- readRDS("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/seu.rds")
modules <- readRDS("project/Liver_Cell_Atlas/Analysis_02/00_Atlas_Markers/modules.rds")


Idents(seu) <- "seurat_clusters"

dir.create("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/signatures") # Create results dir
# Run module score analysis
tictoc::tic("Module Score Analysis")
seu <- module_score_analysis(seu = seu,modules_group = modules, plot_dir = "project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/signatures/", 
                             reduction = "umap", min.cutoff = NA, legend.position = "right", col_pal = NULL, dims_plot = c(1,2), seed = 42,
                             alpha = c(0.1, 0.9), fig.res = 200, pt.size.factor = 1.1, crop = FALSE, raster = FALSE)
tictoc::toc(log = TRUE)


dir.create("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/signatures/ordered_feature") # Create results dir
tictoc::tic("Module Score Print (Ordered Features)")
# Print feature plots with order= TRUE
module_score_print(seu = seu,modules_group = modules, plot_dir = "project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/signatures/ordered_feature/", 
                   reduction = "umap", min.cutoff = NA, legend.position = "right", col_pal = NULL, dims_plot = c(1,2), seed = 42,
                   alpha = c(0.1, 0.9), fig.res = 200, pt.size.factor = 1.1, crop = FALSE, raster = FALSE, order = TRUE)
tictoc::toc(log = TRUE)


saveRDS(seu, "project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/seu.rds")




###############################################################################################

# Here I am going to annotate the cycling cell

seu <- readRDS("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/seu.rds")
Idents(seu) <- "seurat_clusters"

# Read In Lineage Gene Signatures
modules <- readRDS("project/Liver_Cell_Atlas/Analysis_02/00_Atlas_Markers/modules.rds")


cycling <- subset(seu, idents = 18)

covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')

cycling <- cycling %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% ScaleData()
cycling <- RunPCA(cycling, assay="RNA", npcs=40)

set.seed(42)
cycling <- harmony::RunHarmony(cycling, group.by.vars = covariates, reduction = "pca", assay.use = "RNA", reduction.save = "harmony", max.iter.harmony = 50)
cycling <- RunUMAP(cycling, reduction = "harmony", assay = "RNA", dims = 1:20)
cycling <- FindNeighbors(cycling, reduction = "harmony", dims = 1:20)
cycling <- FindClusters(cycling, resolution = 0.5)

dir.create("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/signatures/", recursive = TRUE)

Idents(cycling) <- "seurat_clusters"

png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/Cluster_UMAP.png", width = 8, height = 7, res = 300, units = "in")
print(SeuratPipe::dim_plot(cycling, group.by = "seurat_clusters", label = T, col_pal = SeuratPipe:::discrete_col_pal))
dev.off()


cycling <- module_score_analysis(seu = cycling, modules_group = modules, plot_dir = "project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/signatures/", 
                                 reduction = "umap", min.cutoff = NA, legend.position = "right", col_pal = NULL, dims_plot = c(1,2), seed = 42,
                                 alpha = c(0.1, 0.9), fig.res = 200, pt.size.factor = 1.1, crop = FALSE, raster = FALSE)


saveRDS(cycling, "project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/cycling.rds")


markers <- FindAllMarkers(cycling, only.pos = T, min.pct = 0.25, logfc.threshold = 0.5)
significant_markers <- markers[markers$gene %in% rownames(x = cycling), ] |>
  dplyr::filter(p_val_adj < 0.05) 

CIPR(input_dat = significant_markers,
     comp_method = "logfc_dot_product", 
     reference = "hpca", 
     plot_ind = F,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T)
write.csv(CIPR_top_results, "project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/CIPR_results.csv")

plot <- DotPlot(cycling, features = unique(unlist(modules$Hendo)), group.by = "seurat_clusters") + ggplot2::coord_flip(x=c(1,69), y=c(1,17),clip = "off") + 
  eval(rlang::expr(ggplot2::scale_colour_distiller(name = NULL,type = "div", palette = "RdYlBu"))) + 
  theme(axis.title.y = element_text(margin = margin(r = 120))) +
  
  geom_brace(xstart = 1, xend = 7, ystart = -1.2, yend = -2.5) +
  annotate("text",x=3,y=-3.2,label="Hepatocytes",angle = 45) +
  geom_brace(xstart = 8, xend = 14, ystart = -1.2, yend = -2.5) +
  annotate("text",x=10,y=-3.2,label="Cholangiocytes",angle = 45) +
  geom_brace(xstart = 15, xend = 22, ystart = -1.2, yend = -2.5) +
  annotate("text",x=17.5,y=-3.2,label="Endothelia",angle = 45) +
  geom_brace(xstart = 23, xend = 31, ystart = -1.2, yend = -2.5) +
  annotate("text",x=26,y=-3.2,label="Mesenchyme",angle = 45) +
  geom_brace(xstart = 32, xend = 33, ystart = -1.2, yend = -2.5) +
  annotate("text",x=31.5,y=-3.2,label="Mesothelia",angle = 45) +
  geom_brace(xstart = 34, xend = 37, ystart = -1.2, yend = -2.5) +
  annotate("text",x=34.5,y=-3.2,label="Leucocytes",angle = 45) +
  geom_brace(xstart = 38, xend = 43, ystart = -1.2, yend = -2.5) +
  annotate("text",x=40.5,y=-2.8,label="MPs",angle = 45) +
  geom_brace(xstart = 44, xend = 46, ystart = -1.2, yend = -2.5) +
  annotate("text",x=44,y=-3.2,label="Lymphocytes",angle = 45) +
  geom_brace(xstart = 47, xend = 52, ystart = -1.2, yend = -2.5) +
  annotate("text",x=49,y=-3,label="T Cells",angle = 45) +
  geom_brace(xstart = 53, xend = 54, ystart = -1.2, yend = -2.5) +
  annotate("text",x=53,y=-3,label="B Cells",angle = 45) +
  geom_brace(xstart = 55, xend = 56, ystart = -1.2, yend = -2.5) +
  annotate("text",x=55.5,y=-3,label="ILCs",angle = 45) +
  geom_brace(xstart = 57, xend = 58, ystart = -1.2, yend = -2.5) +
  annotate("text",x=57,y=-3.2,label="Mast Cells",angle = 45) +
  geom_brace(xstart = 59, xend = 60, ystart = -1.2, yend = -2.5) +
  annotate("text",x=59.5,y=-3,label="PDCs",angle = 45) +
  geom_brace(xstart = 61, xend = 64, ystart = -1.2, yend = -2.5) +
  annotate("text",x=61.5,y=-3.2,label="Erythrocytes",angle = 45) +
  geom_brace(xstart = 65, xend = 65, ystart = -1.2, yend = -2.5) +
  annotate("text",x=64,y=-3.2,label="Neutrophils",angle = 45) +
  geom_brace(xstart = 66, xend = 69, ystart = -1.2, yend = -2.5) +
  annotate("text",x=67,y=-3,label="Cycling",angle = 45) +
  
  theme(plot.margin = unit(c(0.01, 0.01, 0.05, 0.05), units="npc"))

png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/Henderson_Module_Gene_by_Cluster_Dotplot.png", width = 15, height = 15, units = "in", res = 300)
print(plot)
dev.off()

# To give annotation to the clusters I am going to write the code line by line incase I want to add note to any annotation
cycling$Alpha.Annotations <- as.character(cycling$seurat_clusters)
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "0"] <- "Hepatocyte"
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "1"] <- "Lymphocyte" # T cell
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "2"] <- "Lymphocyte" # ILC
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "3"] <- "Myeloid Cell" 
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "4"] <- "Lymphocyte"
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "5"] <- "Hepatocyte"
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "6"] <- "Endothelia"
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "7"] <- "Erythrocyte"
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "8"] <- "Lymphocyte" # wasnt as easy to specify due to lower expression of most genes but look to have some lymphocyte expression
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "9"] <- "Myeloid Cell"
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "10"] <- "Lymphocyte" # B cells
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "11"] <- "Myeloid Cell" # cDCs I believe.? Based on CLEC9A expression as well as Atlas Marker expression
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "12"] <- "Myeloid Cell" # Mast Cell
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "13"] <- "Doublet" # has expression of ERG and CD24 as well as Leukocyte markers.
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "14"] <- "Doublet" # has expression of EPCAM and ENG as well as erythrocyte markers, so most likely doublets although a very small amount of expression of EPCAM and ENG is seen in HPA data for erythrocytes.
# Still keeping it as doublets
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "15"] <- "Doublet" # Has hepatocytes, cholangiocyte and mesenchyme genes.
cycling$Alpha.Annotations[cycling$Alpha.Annotations == "16"] <- "Lymphocyte" # B cells
cycling$Alpha.Annotations <- as.factor(cycling$Alpha.Annotations)
table(cycling$Alpha.Annotations)


saveRDS(cycling, "~/project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/cycling.rds")
df <- cycling@meta.data[,"Alpha.Annotations", drop = FALSE]
write.csv(df, "~/project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/cycling_annotations.csv")

# Check how closely the old annotations match the new
tbl <- table(cycling$Alpha.Annotations, cycling$OLD.Alpha.Annotations)
nms <- c("OLD Alpha Annotations","New Alpha Annotations")

# Write the table to a PDF file
pdf("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/Annotation_comparision.pdf", width = 15)
grid::grid.newpage()

tg <- grid::textGrob(nms[1], x=0, hjust=0)
lg <- grid::textGrob(nms[2], x=1, hjust=1)
g <- gridExtra::tableGrob(tbl)
g <- gtable::gtable_add_rows(g, unit(1,"line"), 0)
g <- gtable::gtable_add_cols(g, grid::grobWidth(lg), 0)
g <- gtable::gtable_add_grob(g, tg, t = 1, l=3, r=ncol(g))
g <- gtable::gtable_add_grob(g, lg, l = 1, t=2)
grid::grid.draw(g)
dev.off()
# They seem to match resonable way. Course not 1:1 but nothing stands out massively.

col_pal <- c("Myeloid Cell" = "skyblue2","Mesenchyme" = "darkolivegreen4","Lymphocyte"= "tan3","Hepatocyte" = "#FB9A99","Erythrocyte" = "#CAB2D6", "Endothelia" = "#FDBF6F",
             "Doublets" = "#6699CB","Cholangiocyte" = "indianred", "New Cells" = "darkgrey", "Unknown" = "darkgrey")
png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/Alpha_Annotation_UMAP.png", width = 8, height = 7, res = 300, units = "in")
print(SeuratPipe::dim_plot(cycling, group.by = "Alpha.Annotations", label = T, col_pal = col_pal[levels(cycling$Alpha.Annotations)]))
dev.off()


tmp <- cycling
tmp$seurat_clusters <- factor(tmp$seurat_clusters, levels = c(0,5,   1,2,4,8,10,16,  3,9,11,12, 6, 7,  13,14,15))


plot <- DotPlot(tmp, features = unique(unlist(modules$Hendo)), group.by = "seurat_clusters") + ggplot2::coord_flip(x=c(1,69), y=c(1,17),clip = "off") + 
  eval(rlang::expr(ggplot2::scale_colour_distiller(name = NULL,type = "div", palette = "RdYlBu"))) + 
  theme(axis.title.x = element_text(margin = margin(t = 80))) + theme(axis.title.y = element_text(margin = margin(r = 120))) +
  ggplot2::theme(legend.position = "right") + 
  geom_brace(xstart = -1, xend = -3, ystart = 1, yend = 2,pointing="side") +
  annotate("text",x=-4,y=1.5,label="Hepatocytes") +
  geom_brace(xstart = -1, xend = -3, ystart = 3, yend = 8,pointing="side") +
  annotate("text",x=-4, y= 5.5, label="Lymphocyte") +
  geom_brace(xstart = -1, xend = -3, ystart = 9, yend = 12,pointing="side") +
  annotate("text",x=-4, y = 10.5, label="Myeloid Cell") +
  geom_brace(xstart = -1, xend = -3, ystart = 13, yend = 13,pointing="side") +
  annotate("text",x=-4, y = 13, label="Endothelia") +
  geom_brace(xstart = -1, xend = -4, ystart = 14, yend = 14,pointing="side") +
  annotate("text",x=-5, y = 14, label="Erythrocyte") +
  geom_brace(xstart = -1, xend = -3, ystart = 15, yend = 17,pointing="side") +
  annotate("text",x=-4, y = 16, label="Doublet") +
  
  geom_brace(xstart = 1, xend = 7, ystart = -1.2, yend = -2.5) +
  annotate("text",x=3,y=-3.2,label="Hepatocytes",angle = 45) +
  geom_brace(xstart = 8, xend = 14, ystart = -1.2, yend = -2.5) +
  annotate("text",x=10,y=-3.2,label="Cholangiocytes",angle = 45) +
  geom_brace(xstart = 15, xend = 22, ystart = -1.2, yend = -2.5) +
  annotate("text",x=17.5,y=-3.2,label="Endothelia",angle = 45) +
  geom_brace(xstart = 23, xend = 31, ystart = -1.2, yend = -2.5) +
  annotate("text",x=26,y=-3.2,label="Mesenchyme",angle = 45) +
  geom_brace(xstart = 32, xend = 33, ystart = -1.2, yend = -2.5) +
  annotate("text",x=31.5,y=-3.2,label="Mesothelia",angle = 45) +
  geom_brace(xstart = 34, xend = 37, ystart = -1.2, yend = -2.5) +
  annotate("text",x=34.5,y=-3.2,label="Leucocytes",angle = 45) +
  geom_brace(xstart = 38, xend = 43, ystart = -1.2, yend = -2.5) +
  annotate("text",x=40.5,y=-2.8,label="MPs",angle = 45) +
  geom_brace(xstart = 44, xend = 46, ystart = -1.2, yend = -2.5) +
  annotate("text",x=44,y=-3.2,label="Lymphocytes",angle = 45) +
  geom_brace(xstart = 47, xend = 52, ystart = -1.2, yend = -2.5) +
  annotate("text",x=49,y=-3,label="T Cells",angle = 45) +
  geom_brace(xstart = 53, xend = 54, ystart = -1.2, yend = -2.5) +
  annotate("text",x=53,y=-3,label="B Cells",angle = 45) +
  geom_brace(xstart = 55, xend = 56, ystart = -1.2, yend = -2.5) +
  annotate("text",x=55.5,y=-3,label="ILCs",angle = 45) +
  geom_brace(xstart = 57, xend = 58, ystart = -1.2, yend = -2.5) +
  annotate("text",x=57,y=-3.2,label="Mast Cells",angle = 45) +
  geom_brace(xstart = 59, xend = 60, ystart = -1.2, yend = -2.5) +
  annotate("text",x=59.5,y=-3,label="PDCs",angle = 45) +
  geom_brace(xstart = 61, xend = 64, ystart = -1.2, yend = -2.5) +
  annotate("text",x=61.5,y=-3.2,label="Erythrocytes",angle = 45) +
  geom_brace(xstart = 65, xend = 65, ystart = -1.2, yend = -2.5) +
  annotate("text",x=64,y=-3.2,label="Neutrophils",angle = 45) +
  geom_brace(xstart = 66, xend = 69, ystart = -1.2, yend = -2.5) +
  annotate("text",x=67,y=-3,label="Cycling",angle = 45) +
  
  theme(plot.margin = unit(c(0.01, 0.01, 0.05, 0.05), units="npc"))


png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/Henderson_Module_Gene_by_Cluster_Dotplot_Updated.png", width = 15, height = 15, units = "in", res = 300)
print(plot)
dev.off()



################################################################################################
# Annotations

# Originally I was going to just map previous annotations and any cluster that had 95% of one annotation would be called this. However, looking at the outputs
# There maybe some clusters that will be re defined. I will use these previous annotations to guide most of the new annotations here, but will update some of the clusters.
# Will also provide annotations to all other unknown clusters also and save running a sperate dataset. The only seperate dataset that will need to be run is the cycling cells.


seu <- readRDS("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/seu.rds")

markers <- FindAllMarkers(seu, only.pos = T, min.pct = 0.25, logfc.threshold = 0.5)
significant_markers <- markers[markers$gene %in% rownames(x = seu), ] |>
  dplyr::filter(p_val_adj < 0.05) 

write.csv(markers, "project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/markers.csv")

CIPR(input_dat = significant_markers,
     comp_method = "logfc_dot_product", 
     reference = "hpca", 
     plot_ind = F,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T)

write.csv(CIPR_top_results, "project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/CIPR_results.csv")



# To give annotation to the clusters I am going to write the code line by line incase I want to add note to any annotation
seu$Alpha.Annotations <- as.character(seu$seurat_clusters)
seu$Alpha.Annotations[seu$Alpha.Annotations == "0"] <- "Hepatocyte"
seu$Alpha.Annotations[seu$Alpha.Annotations == "1"] <- "Hepatocyte"
seu$Alpha.Annotations[seu$Alpha.Annotations == "2"] <- "Hepatocyte"
seu$Alpha.Annotations[seu$Alpha.Annotations == "3"] <- "Lymphocyte" # T cell
seu$Alpha.Annotations[seu$Alpha.Annotations == "4"] <- "Myeloid Cell"
seu$Alpha.Annotations[seu$Alpha.Annotations == "5"] <- "Hepatocyte"
seu$Alpha.Annotations[seu$Alpha.Annotations == "6"] <- "Endothelia"
seu$Alpha.Annotations[seu$Alpha.Annotations == "7"] <- "Mesenchyme"
seu$Alpha.Annotations[seu$Alpha.Annotations == "8"] <- "Myeloid Cell" 
seu$Alpha.Annotations[seu$Alpha.Annotations == "9"] <- "Endothelia"
seu$Alpha.Annotations[seu$Alpha.Annotations == "10"] <- "Lymphocyte" # ILC 
seu$Alpha.Annotations[seu$Alpha.Annotations == "11"] <- "Endothelia" 
seu$Alpha.Annotations[seu$Alpha.Annotations == "12"] <- "Lymphocyte" 
seu$Alpha.Annotations[seu$Alpha.Annotations == "13"] <- "Lymphocyte" # Def Lymphocyte - Maybe T cell
seu$Alpha.Annotations[seu$Alpha.Annotations == "14"] <- "Cholangiocyte"
seu$Alpha.Annotations[seu$Alpha.Annotations == "15"] <- "Lymphocyte" # B -Cell
seu$Alpha.Annotations[seu$Alpha.Annotations == "16"] <- "Hepatocyte"
seu$Alpha.Annotations[seu$Alpha.Annotations == "17"] <- "Endothelia"
seu$Alpha.Annotations[seu$Alpha.Annotations == "18"] <- "Cycling"
seu$Alpha.Annotations[seu$Alpha.Annotations == "19"] <- "Lymphocyte" # B cell
seu$Alpha.Annotations[seu$Alpha.Annotations == "20"] <- "Hepatocyte" # potentially myeloid doublet
seu$Alpha.Annotations[seu$Alpha.Annotations == "21"] <- "Endothelia"
seu$Alpha.Annotations[seu$Alpha.Annotations == "22"] <- "Erythrocyte"
seu$Alpha.Annotations[seu$Alpha.Annotations == "23"] <- "Myeloid Cell" # PDCs and cDCs I believe
seu$Alpha.Annotations[seu$Alpha.Annotations == "24"] <- "Mesenchyme"
seu$Alpha.Annotations[seu$Alpha.Annotations == "25"] <- "Hepatocyte"
seu$Alpha.Annotations[seu$Alpha.Annotations == "26"] <- "Hepatocyte"
seu$Alpha.Annotations[seu$Alpha.Annotations == "27"] <- "Hepatocyte"
seu$Alpha.Annotations[seu$Alpha.Annotations == "28"] <- "Myeloid Cell"
seu$Alpha.Annotations[seu$Alpha.Annotations == "29"] <- "Lymphocyte"
seu$Alpha.Annotations[seu$Alpha.Annotations == "30"] <- "Endothelia" # Potentially mesenchyme doublets
seu$Alpha.Annotations[seu$Alpha.Annotations == "31"] <- "Lymphocyte" # B cell
seu$Alpha.Annotations[seu$Alpha.Annotations == "32"] <- "Hepatocyte" # Not convinced it is a Hep but grouped with the Heps but looks more of a doublet. 
                                                                      #Little expression of any lineage markers. Just ASS1 and PECAM1. Also originally was Leukocytes (Alpha.Annotation.Tmp) but changed to Heps
seu$Alpha.Annotations[seu$Alpha.Annotations == "33"] <- "Myeloid Cell" # Most likely a myeloid populations according to markers from the HCA list. # Previously was endothelia in Alpha.Annotation.Tmp
seu$Alpha.Annotations[seu$Alpha.Annotations == "34"] <- "Hepatocyte" # Potential doublet
seu$Alpha.Annotations[seu$Alpha.Annotations == "35"] <- "Hepatocyte" # Most likely erythrocyte doublets
seu$Alpha.Annotations[seu$Alpha.Annotations == "36"] <- "Hepatocyte" # Change from Myeloid Cells (Alpha.Annotation.Tmp). Potentially B cell doublets.
seu$Alpha.Annotations <- as.factor(seu$Alpha.Annotations)
table(seu$Alpha.Annotations)




# Check how closely the old annotations match the new
tbl <- table(seu$Alpha.Annotations, seu$OLD.Alpha.Annotations)
nms <- c("OLD Alpha Annotations","New Alpha Annotations")

# Write the table to a PDF file
pdf("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Annotation_comparision.pdf", width = 15)
grid::grid.newpage()

tg <- grid::textGrob(nms[1], x=0, hjust=0)
lg <- grid::textGrob(nms[2], x=1, hjust=1)
g <- gridExtra::tableGrob(tbl)
g <- gtable::gtable_add_rows(g, grid::unit(1,"line"), 0)
g <- gtable::gtable_add_cols(g, grid::grobWidth(lg), 0)
g <- gtable::gtable_add_grob(g, tg, t = 1, l=3, r=ncol(g))
g <- gtable::gtable_add_grob(g, lg, l = 1, t=2)
grid::grid.draw(g)
dev.off()
# They seem to match resonable way. Course not 1:1 but nothing stands out massively.

col_pal <- c("Myeloid Cell" = "skyblue2","Mesenchyme" = "darkolivegreen4","Lymphocyte"= "tan3","Hepatocyte" = "#FB9A99","Erythrocyte" = "#CAB2D6", "Endothelia" = "#FDBF6F",
             "Doublets" = "#6699CB","Cholangiocyte" = "indianred", "New Cells" = "darkgrey", "Unknown" = "darkgrey")
png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Alpha_Annotation_UMAP.png", width = 8, height = 7, res = 300, units = "in")
print(SeuratPipe::dim_plot(seu, group.by = "Alpha.Annotations", label = T, col_pal = col_pal[levels(seu$Alpha.Annotations)], raster = FALSE))
dev.off()

# Read In Lineage Gene Signatures
modules <- readRDS("project/Liver_Cell_Atlas/Analysis_02/00_Atlas_Markers/modules.rds")


# Re create dot plot with genes with new Alpha Annotations
tmp <- seu
tmp$seurat_clusters <- factor(seu$seurat_clusters, levels = c(0,1,2,5,16,20,25,26,27,32,34,35,36,3,10,12,13,15,19,29,31,4,8,23,28,33,6,9,11,17,21,30,7,24,14,22,18))

plot <- DotPlot(tmp, features = unique(unlist(modules$Hendo)), group.by = "seurat_clusters") + ggplot2::coord_flip(x=c(1,69), y=c(1,37),clip = "off") + 
  eval(rlang::expr(ggplot2::scale_colour_distiller(name = NULL,type = "div", palette = "RdYlBu"))) + 
  theme(axis.title.x = element_text(margin = margin(t = 80))) + theme(axis.title.y = element_text(margin = margin(r = 120))) +
  ggplot2::theme(legend.position = "right") + 
  geom_brace(xstart = -1, xend = -3, ystart = 1, yend = 13,pointing="side") +
  annotate("text",x=-5,y=6,label="Hepatocytes",angle = 45) +
  geom_brace(xstart = -1, xend = -3, ystart = 14, yend = 21,pointing="side") +
  annotate("text",x=-5,y=16.5,label="Lymphocyte",angle = 45) +
  geom_brace(xstart = -1, xend = -3, ystart = 22, yend = 26,pointing="side") +
  annotate("text",x=-5,y=23,label="Myeloid Cell",angle = 45) +
  geom_brace(xstart = -1, xend = -3, ystart = 27, yend = 32,pointing="side") +
  annotate("text",x=-4.8,y=28.7,label="Endothelia",angle = 45) +
  geom_brace(xstart = -1, xend = -3, ystart = 33, yend = 34,pointing="side") +
  annotate("text",x=-5,y=32.2,label="Mesenchyme",angle = 45) +
  geom_brace(xstart = -1, xend = -3, ystart = 35, yend = 35,pointing="side") +
  annotate("text",x=-5.1,y=33.8,label="Cholangiocyte",angle = 45) +
  geom_brace(xstart = -1, xend = -3, ystart = 36, yend = 36,pointing="side") +
  annotate("text",x=-4.9,y=35,label="Erythrocyte",angle = 45) +
  geom_brace(xstart = -1, xend = -3, ystart = 37, yend = 37,pointing="side") +
  annotate("text",x=-4.5,y=36.3,label="Cycling",angle = 45) +
  
  geom_brace(xstart = 1, xend = 7, ystart = -3, yend = -5) +
  annotate("text",x=3,y=-6.5,label="Hepatocytes",angle = 45) +
  geom_brace(xstart = 8, xend = 14, ystart = -3, yend = -5) +
  annotate("text",x=10.5,y=-6.5,label="Cholangiocytes",angle = 45) +
  geom_brace(xstart = 15, xend = 22, ystart = -3, yend = -5) +
  annotate("text",x=17.5,y=-6.5,label="Endothelia",angle = 45) +
  geom_brace(xstart = 23, xend = 31, ystart = -3, yend = -5) +
  annotate("text",x=26,y=-6.5,label="Mesenchyme",angle = 45) +
  geom_brace(xstart = 32, xend = 33, ystart = -3, yend = -5) +
  annotate("text",x=31.5,y=-6.5,label="Mesothelia",angle = 45) +
  geom_brace(xstart = 34, xend = 37, ystart = -3, yend = -5) +
  annotate("text",x=34.5,y=-6.5,label="Leucocytes",angle = 45) +
  geom_brace(xstart = 38, xend = 43, ystart = -3, yend = -5) +
  annotate("text",x=40.5,y=-5.7,label="MPs",angle = 45) +
  geom_brace(xstart = 44, xend = 46, ystart = -3, yend = -5) +
  annotate("text",x=44,y=-6.5,label="Lymphocytes",angle = 45) +
  geom_brace(xstart = 47, xend = 52, ystart = -3, yend = -5) +
  annotate("text",x=49,y=-6.2,label="T Cells",angle = 45) +
  geom_brace(xstart = 53, xend = 54, ystart = -3, yend = -5) +
  annotate("text",x=53,y=-6.2,label="B Cells",angle = 45) +
  geom_brace(xstart = 55, xend = 56, ystart = -3, yend = -5) +
  annotate("text",x=55.5,y=-5.7,label="ILCs",angle = 45) +
  geom_brace(xstart = 57, xend = 58, ystart = -3, yend = -5) +
  annotate("text",x=57,y=-6.3,label="Mast Cells",angle = 45) +
  geom_brace(xstart = 59, xend = 60, ystart = -3, yend = -5) +
  annotate("text",x=59.5,y=-6,label="PDCs",angle = 45) +
  geom_brace(xstart = 61, xend = 64, ystart = -3, yend = -5) +
  annotate("text",x=61.5,y=-6.5,label="Erythrocytes",angle = 45) +
  geom_brace(xstart = 65, xend = 65, ystart = -3, yend = -5) +
  annotate("text",x=64,y=-6.5,label="Neutrophils",angle = 45) +
  geom_brace(xstart = 66, xend = 69, ystart = -3, yend = -5) +
  annotate("text",x=67,y=-6,label="Cycling",angle = 45) +
  
  theme(plot.margin = unit(c(0.01, 0.01, 0.05, 0.05), units="npc"))

png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Henderson_Module_Gene_by_Cluster_Dotplot_Updated.png", width = 15, height = 15, units = "in", res = 300)
print(plot)
dev.off()




# Read in cycling cell annotation csv and map to the dataset
cycling_annoations <- read.csv("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/cycling/cycling_annotations.csv", row.names = 1)


names <- colnames(seu)[seu$Alpha.Annotations == "Cycling"]
seu$Alpha.Annotations <- as.character(seu$Alpha.Annotations)
seu$Alpha.Annotations[rownames(cycling_annoations)] <- cycling_annoations[rownames(cycling_annoations),"Alpha.Annotations"]
#Check they have mapped ok
df <- seu@meta.data[names, "Alpha.Annotations", drop = FALSE]
all.equal(df, cycling_annoations)

# make alpha annotation a factor
seu$Alpha.Annotations <- as.factor(seu$Alpha.Annotations)

col_pal <- c("Myeloid Cell" = "skyblue2","Mesenchyme" = "darkolivegreen4","Lymphocyte"= "tan3","Hepatocyte" = "#FB9A99","Erythrocyte" = "#CAB2D6", "Endothelia" = "#FDBF6F",
             "Doublets" = "#6699CB","Cholangiocyte" = "indianred", "New Cells" = "darkgrey", "Unknown" = "darkgrey")
png("project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Alpha_Annotation_UMAP_Annotated_Cycling.png", width = 8, height = 7, res = 300, units = "in")
print(SeuratPipe::dim_plot(seu, group.by = "Alpha.Annotations", label = T, col_pal = col_pal[levels(seu$Alpha.Annotations)], raster = FALSE))
dev.off()


saveRDS(seu, "~/project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/seu.rds")
df <- seu@meta.data[,"Alpha.Annotations", drop = FALSE]
write.csv(df, "~/project/Liver_Cell_Atlas/Analysis_02/01_Broad_Cell_Annotations/results/Alpha_Annotations.csv")
write.csv(df, "~/project/Liver_Cell_Atlas/Analysis_02/01_Alpha_Annotations.csv")
