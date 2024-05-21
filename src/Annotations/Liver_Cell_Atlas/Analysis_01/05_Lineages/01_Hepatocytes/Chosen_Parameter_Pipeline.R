# Name



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
seu <- readRDS("~/project/Liver_Cell_Atlas/Analysis_01/01_Broad_Cell_Annotations/seu.rds") 
opts <- list()
opts$modules_group <- readRDS("~/project/Liver_Cell_Atlas/Analysis_01/00_Atlas_Markers/modules.rds")



## Subset Lineage

##-----------------------------------------------------------------------
# Subset to retain cells from specified lineage
##-----------------------------------------------------------------------
opts$lineage_subset <- c("Hepatocyte") #subset lineage

seu_anno <- read.csv(
  file = paste0("~/project/Liver_Cell_Atlas/Analysis_01/Annotations_4.csv")) |> # Read in cluster information
  tibble::column_to_rownames(var = "X") |>
  dplyr::select(c(Preliminary.Annotations, Alpha.Annotations, Potential.Doublets)) |> 
  dplyr::mutate_at(c("Preliminary.Annotations", "Alpha.Annotations", "Potential.Doublets"), factor)
# Add metadata information to current Seurat object and filter
seu <- AddMetaData(seu, metadata = seu_anno)
Idents(object = seu) <- "Alpha.Annotations"
seu <- subset(seu, idents = c(opts$lineage_subset), invert = FALSE)

dim_plot(seu, group.by = "Alpha.Annotations", label = T)






## Clean Seurat Object
DefaultAssay(seu) <- "RNA"

#Trim down seurat object using diet seurat
seu <- DietSeurat(seu, assays = "RNA")
#Remove metadata associated with previous clustering and annotations
seu@meta.data <- seu@meta.data[, -which(colnames(seu@meta.data) %in% c(grep(colnames(seu@meta.data), pattern = "_snn_", value = T),grep(colnames(seu@meta.data), pattern = "^seurat_clusters", value = T)))]
#Remove previous signature scores
seu@meta.data <- seu@meta.data[, -which(colnames(seu@meta.data) %in% gsub("^.*?\\.", "", names(flatten(opts$modules_group))))]
#Remove variable genes
VariableFeatures(seu) <- NULL


## Fix Metadata

# Fixed the NAs in Ethnicity
seu$donor_ethnicity[is.na(seu$donor_ethnicity)] <- "unknown"




# Print MetaData
Hmisc::describe(seu@meta.data)




#  SeuratPipe

# Opts

# Path to folder containing this script
path <- "project/Liver_Cell_Atlas/Analysis_01/05_Lineages/01_Hepatocytes/"
# Name you want the output results folder
results_folder <- "chosen_parameter" # eg chosen_parameters
# Name you want to the sinked console output txt file
sink_name <- "chosen_parameter" # console_output



# Add in sample as metadata with the correct case for the pipeline
seu$sample <- as.factor(seu$SAMPLE)
# Workflow needs condition variable to be present
seu$condition <- as.factor("Healthy")


# Number of PCs, can be a vector: c(30, 50)
opts$npcs <- c(40)
# Top PCs to perform UMAP and clustering, can be a vector: c(30, 40)
opts$ndims <- c(40)
# Clustering resolutions
opts$res <- seq(from = 0.5, to = 0.5, by = 0.1)
# Batch ID to perform integration
opts$batch_id <- c("study", "SAMPLE")
# QC information to plot
opts$qc_to_plot <- c("nCount_RNA","nFeature_RNA", "percent.mt", "remain_counts_soupX","pct_contamin")
# Metadata columns to plot
opts$meta_to_plot <- c("suspension_type", "scmap_cluster_anno","scmap_cell_anno", "Phase","general_labs","consistent_labs","marker_labs","marker_general_labs","assay",
                       "tissue","donor_sex","donor_ethnicity","study","customclassif","Preliminary.Annotations", "Alpha.Annotations", "Potential.Doublets")


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
opts$n_hvgs <- 3000
# Set specific seed for reproducibility
opts$seed <- 42
# Discrete colour palette
opts$discrete_col_pal <- SeuratPipe:::discrete_col_pal
# Whether to label the clusters in 'plot_reduction' space.
opts$label <- TRUE
# Sets size of labels.
opts$label.size <- 6
# Adjust point size for plotting.
opts$pt.size <- 1.4
# Figure resolution in ppi
opts$fig.res = 200


# Load modules
# see above

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


sink()
sink(type = "message")


# Add Automated Annoations

# Make Grouped Marker General Annotation
Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}


seu$marker_labs_group <- NA
for(i in levels(seu$seurat_clusters)){
  txt <- Modes(seu$marker_labs[seu$seurat_clusters==i])
  seu$marker_labs_group[seu$seurat_clusters==i] <- txt
}

# ScType Re run
lapply(c("HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# Use marker gene list from liver group 
marker_genes = read.csv("~/project/Liver_Cell_Atlas/Analysis_01/00_Atlas_Markers/marker_lists/Atlas_Marker_List_Collapsed_Updated.csv", header = T, sep = ',')
marker_list = split(marker_genes$Gene, marker_genes$Specific.Type)
names(marker_list)
gs_neg <- vector("list", length(marker_list));names(gs_neg) <- names(marker_list)
es.max = sctype_score(scRNAseqData = seu[["RNA"]]@scale.data, scaled=TRUE, gs=marker_list, gs2=gs_neg)


# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seu@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seu@meta.data[seu@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

seu@meta.data$customclassif_new = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seu@meta.data$customclassif_new[seu@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

dim_plot(seu, group.by = "customclassif_new", raster = FALSE)

# Save Data

# Save updated metadata
opts$preclean_annoation_dir <- io$out_dir
saveRDS(object = list(seu = seu, opts = opts),
        file = paste0(opts$preclean_annoation_dir,"/", "Hepatocytes_preclean.rds"))

file.remove(paste0(io$out_dir,"/",grep(paste0(opts$obj_filename,"*.*rds"), list.files(io$out_dir), value = TRUE))) # Remove intermediate RDS files







# Presenation Documents
#Make sure document is saved before running this
# Pull Results dir from eariler in script
this_doc <- "project/Liver_Cell_Atlas/Analysis_01/05_Lineages/01_Hepatocytes/Chosen_Parameter_Pipeline.R" #fill in name of this document
lines <- grep("*Input/Output*",readLines(paste0(this_doc)))
text <- readLines(this_doc)[c(lines[1],lines[1]+1)]
eval(parse(text=text))

# Produce document (needs io$our_dir, code above fetches this from earlier)
presentation_create(title = "Liver Cell Atlas - Hepatocytes - Data Presentation",output_file = "~/project/Liver_Cell_Atlas/Analysis_01/05_Lineages/01_Hepatocytes/Hepatocyte_Data_Presenation.html")



