# Liver Cell Atlas - Cholangiocyte (Clean) - Chosen Parameters
  


# Load Libraries and Environment Setup - Data Preparation
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratPipe))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
source("~/project/Functions/Utils.R")
seed <- 42 # Seed number for scipt
set.seed(seed) 


# Data Preparation
## Load Data

obj <- readRDS("~/project/Liver_Cell_Atlas/Analysis_02/02_Lineages_itr01/03_Cholangiocytes/01_PreClean/chosen_parameters/Cholangiocyte_unclean.rds") 
seu <- obj$seu
opts <- list()
opts$modules_group <- readRDS("~/project/Liver_Cell_Atlas/Analysis_02/00_Atlas_Markers/modules.rds")



## Subset Lineage

##-----------------------------------------------------------------------
# Subset to retain cells from specified lineage
##-----------------------------------------------------------------------
opts$lineage_anno <- "Cholangiocyte" #subset lineage

Idents(object = seu) <- "Beta.Annotation.Lineage"
seu <- subset(seu, idents = c(opts$lineage_anno), invert = FALSE)

dim_plot(seu, group.by = "Beta.Annotation.Lineage", label = T)






## Clean Seurat Object
DefaultAssay(seu) <- "RNA"
#Trim down seurat object using diet seurat
seu <- DietSeurat(seu, assays = "RNA")
#Remove metadata associated with previous clustering and annotations
seu@meta.data <- seu@meta.data[, -which(colnames(seu@meta.data) %in% c(grep(colnames(seu@meta.data), pattern = "_snn_", value = T),grep(colnames(seu@meta.data), pattern = "^seurat_clusters", value = T)))]
#Remove previous signature scores
seu@meta.data <- seu@meta.data[, -which(colnames(seu@meta.data) %in% gsub("^.*?\\.", "", names(flatten(opts$modules_group))))]
#Remove Previous Cluster scores.
seu@meta.data <- seu@meta.data[, -which(colnames(seu@meta.data) %in%  c(grep(colnames(seu@meta.data), pattern = "^Cluster", value = T)))]
# Remove Malat Threshold, and sc-type annotations from previous analysis
seu$malat1_threshold <- NULL
seu$HCA_General_sctype <- NULL
seu$HCA_Specific_sctype <- NULL
seu$Hendo_sctype <- NULL
#Remove variable genes
VariableFeatures(seu) <- NULL


# Print MetaData
Hmisc::describe(seu@meta.data)




#  SeuratPipe

# Opts

# Path to folder containing this script
path <- "~/project/Liver_Cell_Atlas/Analysis_02/02_Lineages_itr01/03_Cholangiocytes/02_Clean/"
# Name you want the output results folder
results_folder <- "chosen_parameters" # eg chosen_parameters
# Name you want to the sinked console output txt file
sink_name <- "02_Cholangiocyte_Clean_Chosen_Parameters_console" # console_output
covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')


# Add in sample as metadata, this isnt the same as SAMPLE but is the same as study. The reason for this is the need to a large amount of colours for the pipeline to work otherwise. Plots looking at SAMPLE can be made afterwards.
seu$sample <- as.factor(seu$STUDY)
# Workflow needs condition variable to be present
seu$condition <- as.factor("Healthy")


# Number of PCs, can be a vector: c(30, 50)
opts$npcs <- c(20)
# Top PCs to perform UMAP and clustering, can be a vector: c(30, 40)
opts$ndims <- c(20)
# Clustering resolutions
opts$res <- seq(from = 0.2, to = 0.2, by = 0.1)
# Batch ID to perform integration
opts$batch_id <- covariates
# QC information to plot
opts$qc_to_plot <- c("nCount_RNA","nFeature_RNA", "percent.mt")
# Metadata columns to plot
opts$meta_to_plot <-  c("suspension_type", "Phase","assay","donor_age",
                        "tissue","donor_sex","donor_ethnicity","STUDY","Potential.Doublets")


# Test genes that are detected in a minimum fraction of min.pct cells
opts$min.pct <- 0.25
# Test genes that show, on average, at least X-fold difference
# between the two groups of cells.
opts$logfc.threshold <- 0.5
# Only return positive markers
opts$only.pos <- TRUE
# Maximum number of marker genes to plot for each cluster
opts$topn_genes <- 10
# Retain marker genes per cluster if their
# `pct.1 - pct.2 > diff_cluster_pct`, i.e. they show cluster
# specific expression. Set to -Inf, to ignore this additional filtering.
opts$diff_cluster_pct <- 0.1
# Adjusted p-value threshold to consider marker genes per cluster.
opts$pval_adj <- 0.05
# Should we create feature plots of cluster markers?
opts$plot_cluster_markers = TRUE
# Which PCs to remove prior to running harmony (technical effects)
# Can be a vector e.g. c(1, 2, 4). If NULL, all PCs are used as input.
opts$pcs_to_remove <- NULL
# Maximum cutoff values for plotting continuous features, e.g. gene expression
# Gives better plots where colour scale is not driven by a (few) outlier cells.
# Set to NULL to recoved default Seurat plots
opts$max.cutoff <- "q98"
# Filename of the integrated object to be stored on disk
opts$obj_filename <- "seu_harmony"
# If integrated Harmony file exists, should we force re-analysis
# of Harmony, or read object? For computing time efficiency purposes.
opts$force_reanalysis = TRUE


# Number of highly variable genes to compute
opts$n_hvgs <- 2000
# Set specific seed for reproducibility
opts$seed <- 42
# Discrete colour palette
opts$discrete_col_pal <- c("#D26263", "#6399CA", "#FABD6A", "#CBB3D8", "#FD9D9C", "#CD8240", "#6D8938", "#A7A6A7", "#7CBEEB", "#A883FD", "#9ACB9A", "#EEE484", 
                           "#CBCBC0", "#2A6089", "#683D97", "#329460", "#FF85FC", "#2226F6", "#FE1692", "#FED626", "#1CCFD2", "#8F451C", "#FE7E00", "#228EFD",
                           "#D0CF0D", "#DE65FE", "#942622", "#90846D", "#B53568", "#5C565C", "#FD222A", "#EEC790", "#26FBFB", "#BFF26C", "#9BD035", "#0DFC35",
                           "#00FFB2", "#FB00CB", "#890076", "#FB9CD4", "#BA001C", "#A500FA", "#A36375", "#FE0DFA", "#D39EF5", "#FE60BA", "#22817E" )
# Whether to label the clusters in 'plot_reduction' space.
opts$label <- TRUE
# Sets size of labels.
opts$label.size <- 6
# Adjust point size for plotting.
opts$pt.size <- 1.4
# Figure resolution in ppi
opts$fig.res = 200



# Load modules
# See above



##-------------------------------------------------------
# Check Colour Numbers
if (is.list(seu)){
  for (i in seu){
    if (max(sapply(lapply(i@meta.data[,opts$meta_to_plot], factor), nlevels)) > length(opts$discrete_col_pal)){
      stop("Need more colours in discrete_col_pal for the number of levels in the selected metadata")
    }
  }
} else {
  if (max(sapply(lapply(seu@meta.data[,opts$meta_to_plot], factor), nlevels)) > length(opts$discrete_col_pal)){
    stop("Need more colours in discrete_col_pal for the number of levels in the selected metadata")
  }
}
##--------------------------------------------------------


##--------------------------------------------------------
# Makes user press enter to continue when PCs are removed to check the right number of PCs are given for ndims to account for missing PCs
readkey <- function() {
  line <- readline(prompt="Press [enter] to continue")
}
if (!is.null(opts$pcs_to_remove)){
  print("Stop Loop if you haven't changed the numbers in ndims to account for any # of missing PCs")
  readkey()
}
##-----------------------------------------------------

# Output folder
io <- list() # Input/Output
io$out_dir <- paste0(path,"/",results_folder) # Name of result folder 
if (!dir.exists(io$out_dir)) dir.create(io$out_dir, recursive = TRUE)

# Sink console output

file <- file(paste0(path, "/", sink_name, ".txt"))
sink(file = file, append = TRUE, split = TRUE)
sink(file = file, append = TRUE, type = "message")
Sys.Date()
Sys.time()
print("") # New line
print(paste0("Outs folder: ",path,"/",results_folder)) # print outs folder to sink
print("") # new line
print(opts) # print opts to sink
print("") # new line
##-----------------------------------------------------------------------
# Analysis pipeline
##-----------------------------------------------------------------------

# Perform data integration using Harmony
seu <- run_harmony_pipeline(
  seu_obj = seu,
  out_dir = io$out_dir,
  batch_id = opts$batch_id,
  npcs = opts$npcs,
  ndims = opts$ndims,
  res = opts$res,
  modules_group = opts$modules_group,
  metadata_to_plot = opts$meta_to_plot,
  qc_to_plot = opts$qc_to_plot,
  logfc.threshold = opts$logfc.threshold,
  min.pct = opts$min.pct,
  only.pos = opts$only.pos,
  topn_genes = opts$topn_genes,
  diff_cluster_pct = opts$diff_cluster_pct,
  pval_adj = opts$pval_adj,
  pcs_to_remove = opts$pcs_to_remove,
  obj_filename = opts$obj_filename,
  force_reanalysis = opts$force_reanalysis,
  plot_cluster_markers = opts$plot_cluster_markers,
  max.cutoff = opts$max.cutoff,
  n_hvgs = opts$n_hvgs,
  seed = opts$seed,
  discrete_col_pal = opts$discrete_col_pal,
  label = opts$label,
  label.size = opts$label.size,
  pt.size = opts$pt.size,
  fig.res = opts$fig.res,
  heatmap_downsample_cols = 500,
  raster = FALSE)

sessionInfo()


sink()
sink(type = "message")




# Add Automated Annoations


# Load ScType and make tmp data for scaling all genes 
lapply(c("HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
tmp <- seu
tmp <- ScaleData(tmp, features = rownames(tmp))

#Run ScType On HCA General markers
# Use marker gene list from liver group 
marker_list <- opts$modules_group$Atlas_general
gs_neg <- vector("list", length(marker_list));names(gs_neg) <- names(marker_list)
es.max = sctype_score(scRNAseqData = tmp[["RNA"]]@scale.data, scaled=TRUE, gs=marker_list, gs2=gs_neg)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(tmp@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(tmp@meta.data[tmp@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(tmp@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

seu@meta.data$HCA_General_sctype = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seu@meta.data$HCA_General_sctype[seu@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


#Run ScType On HCA Specific markers
# Use marker gene list from liver group 
marker_list <- opts$modules_group$Atlas_specific
gs_neg <- vector("list", length(marker_list));names(gs_neg) <- names(marker_list)
es.max = sctype_score(scRNAseqData = tmp[["RNA"]]@scale.data, scaled=TRUE, gs=marker_list, gs2=gs_neg)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(tmp@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(tmp@meta.data[tmp@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(tmp@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

seu@meta.data$HCA_Specific_sctype = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seu@meta.data$HCA_Specific_sctype[seu@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


#Run ScType On HCA Specific markers
# Use marker gene list from liver group 
marker_list <- opts$modules_group$Hendo
gs_neg <- vector("list", length(marker_list));names(gs_neg) <- names(marker_list)
es.max = sctype_score(scRNAseqData = tmp[["RNA"]]@scale.data, scaled=TRUE, gs=marker_list, gs2=gs_neg)

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(tmp@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(tmp@meta.data[tmp@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(tmp@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

seu@meta.data$Hendo_sctype = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seu@meta.data$Hendo_sctype[seu@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


dim_plot(seu, group.by = "HCA_General_sctype")
dim_plot(seu, group.by = "HCA_Specific_sctype")
dim_plot(seu, group.by = "Hendo_sctype")



# Malat Thresholding
# This script filters cells based on MALAT1 expression
# Code is written assuming the use of a Seurat object, but
# should be able to be applied to any single-cell object

# Look at histogram of normalized MALAT1 expression
hist(seu@assays$RNA@data["MALAT1",], freq = FALSE, breaks = 100)


# Run this function on MALAT1 reads, eg:
threshold <- define_malat1_threshold(seu@assays$RNA@data["MALAT1",], max_counts = 2)
opts$malat1_threshold <- threshold

# OR assign TRUE/FALSE values to determine which cells would pass the filter
malat1_threshold <- seu@assays$RNA@data['MALAT1',] > threshold
seu$malat1_threshold <- malat1_threshold
DimPlot(seu, group.by = "malat1_threshold")



annotations <- list(Beta.Annotation.Lineage = c("0" = "Cholangiocyte",
                                                "1" = "Cholangiocyte",
                                                "2" = "Cholangiocyte",
                                                "3" = "Cholangiocyte",
                                                "4" = "Cholangiocyte",
                                                "5" = "Cholangiocyte",
                                                "6" = "Cholangiocyte",
                                                "7" = "Cholangiocyte",
                                                "8" = "Doublet"),
                    Beta.Annotation.SubLineage = c("0" = "ApoLipo / Small",
                                                   "1" = "Keratin / Small / Ribosomal / Cycling",
                                                   "2" = "Mucus Secreting / Large / Keratin",
                                                   "3" = "ApoLipo / Small",
                                                   "4" = "LAMC2+ / Small",
                                                   "5" = "Stressed?",
                                                   "6" = "Keratin / CXCL8+",
                                                   "7" = "Small (Keratin)",
                                                   "8" = "TBC"),
                    Potential.Doublets.2 = c("0" = "TRUE",
                                             "1" = "FALSE",
                                             "2" = "TRUE",
                                             "3" = "TRUE",
                                             "4" = "FALSE",
                                             "5" = "FALSE",
                                             "6" = "FALSE",
                                             "7" = "FALSE",
                                             "8" = "TRUE"))






seu$Beta.Annotation.Lineage <- plyr::mapvalues(seu$seurat_clusters, from = levels(seu$seurat_clusters), annotations$Beta.Annotation.Lineage)
table(seu$Beta.Annotation.Lineage)

seu$Beta.Annotation.SubLineage <- plyr::mapvalues(seu$seurat_clusters, from = levels(seu$seurat_clusters), annotations$Beta.Annotation.SubLineage)
table(seu$Beta.Annotation.SubLineage)

seu$Potential.Doublets.2 <- plyr::mapvalues(seu$seurat_clusters, from = levels(seu$seurat_clusters), annotations$Potential.Doublets.2)
table(seu$Potential.Doublets.2)

DimPlot(seu, group.by = "Beta.Annotation.Lineage")
DimPlot(seu, group.by = "Beta.Annotation.SubLineage")
DimPlot(seu, group.by = "Potential.Doublets.2")

df <- seu@meta.data[,c("Beta.Annotation.Lineage","Beta.Annotation.SubLineage","Potential.Doublets","Potential.Doublets.2" )]
write.csv(df, "~/project/Liver_Cell_Atlas/Analysis_02/02_Lineages_itr01/03_Cholangiocytes/02_Clean/chosen_parameters/Cholangiocyte_Clean_annotations.csv")

cholangiocyte_df <- read.csv("~/project/Liver_Cell_Atlas/Analysis_02/02_Lineages_itr01/03_Cholangiocytes/01_PreClean/chosen_parameters/Cholangiocyte_unclean_annotations.csv", row.names = 1)
cholangiocyte_df[names(seu$Beta.Annotation.Lineage),"Beta.Annotation.Lineage"] <- as.character(seu$Beta.Annotation.Lineage)
cholangiocyte_df$Beta.Annotation.SubLineage <- "TBC"
cholangiocyte_df[names(seu$Beta.Annotation.SubLineage),"Beta.Annotation.SubLineage"] <- as.character(seu$Beta.Annotation.SubLineage)
cholangiocyte_df$Potential.Doublets.2 <- "FALSE"
cholangiocyte_df[names(seu$Potential.Doublets.2),"Potential.Doublets.2"] <- as.character(seu$Potential.Doublets.2)


# Test that the above re naming is correct
test <- seu
test$Beta.Annotation.Lineage <- NULL
test$Potential.Doublets.2 <- NULL
test$Beta.Annotation.SubLineage <- NULL
test <- AddMetaData(test, cholangiocyte_df)

all.equal(as.character(test$Beta.Annotation.Lineage), as.character(seu$Beta.Annotation.Lineage))
all.equal(as.character(test$Beta.Annotation.SubLineage), as.character(seu$Beta.Annotation.SubLineage))
all.equal(as.character(test$Potential.Doublets.2), as.character(seu$Potential.Doublets.2))
DimPlot(test, group.by = "Beta.Annotation.Lineage")
DimPlot(test, group.by = "Beta.Annotation.SubLineage")
DimPlot(test, group.by = "Potential.Doublets.2")

write.csv(cholangiocyte_df, "~/project/Liver_Cell_Atlas/Analysis_02/02_Lineages_itr01/03_Cholangiocytes/Cholangiocyte_annotations.csv")

# Save Data
# Save updated metadata
opts$Cholangiocyte_clean_anno_dir <- io$out_dir
saveRDS(object = list(seu = seu, opts = opts),
        file = paste0(opts$Cholangiocyte_clean_anno_dir,"/", "Cholangiocyte_clean.rds"))

file.remove(paste0(io$out_dir,"/",grep(paste0(opts$obj_filename,"*.*rds"), list.files(io$out_dir), value = TRUE))) # Remove intermediate RDS files


# Presenation Documents
# Produce document (needs io$our_dir, code above fetches this from earlier)
presentation_create(title = "Liver Cell Atlas - Cholangiocyte (Clean ish) - Data Presentation", 
                    template = "~/project/Liver_Cell_Atlas/Analysis_02/02_Lineages_itr01/03_Cholangiocytes/02_Clean/Cholangiocyte_Clean_template.Rmd", 
                    obj_path = "~/project/Liver_Cell_Atlas/Analysis_02/02_Lineages_itr01/03_Cholangiocytes/02_Clean/chosen_parameters/Cholangiocyte_clean.rds", 
                    results_path = "~/project/Liver_Cell_Atlas/Analysis_02/02_Lineages_itr01/03_Cholangiocytes/02_Clean/chosen_parameters/",
                    final_cluster_param = c("Beta.Annotation.Lineage","Beta.Annotation.SubLineage", "Potential.Doublets.2"),
                    output_file = "~/project/Liver_Cell_Atlas/Analysis_02/02_Lineages_itr01/03_Cholangiocytes/02_Clean/Cholangiocyte_Clean_Data_Presenation.html")



### Save Harmony Components
obj <- readRDS("~/project/Liver_Cell_Atlas/Analysis_02/02_Lineages_itr01/03_Cholangiocytes/02_Clean/chosen_parameters/Cholangiocyte_clean.rds")
seu <- obj$seu
opts <- obj$opts
path <- paste0(opts$Cholangiocyte_clean_anno_dir, "/harmony_backup/")

dir.create(path)

write.csv(seu@reductions$harmony@cell.embeddings, paste0(path, "cell_embeddings.csv"))
write.csv(seu@reductions$harmony@feature.loadings, paste0(path, "feature_loadings.csv"))
write.csv(seu@reductions$harmony@feature.loadings.projected, paste0(path, "feature_loadings_projected.csv"))


