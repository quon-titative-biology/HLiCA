# Basic clustering with Seurat
library(Seurat)
library(dplyr)
library(scClustViz)
library(viridis)
library(presto)
library(cluster)
library(viridisLite)
library(shiny)
library(DropletQC)
library(dplyr)
library(ggplot2)

woodchuck <- "3391"
sample <- "PBMC"

# Get woodchuck metadata
source("~/Dropbox/Zoe/scf_version/analysis/scripts/woodchuckMetadata.R")

# Isolate metadata for specific sample
sampleInfo <- woodchuck_info[woodchuck_info$woodchuck == woodchuck & 
                               woodchuck_info$sample_type == sample,]

# Load in data
sobj.data <- Read10X(paste("~/temp_cellranger_outputs/",
                           woodchuck, "/", sample, "/filtered_feature_bc_matrix/", sep = ""))
sobj <- CreateSeuratObject(counts = sobj.data, min.features = 200, min.cells = 3)

if(sample != "PBMC") {
  sobj <- AddMetaData(sobj, sampleInfo$orig.ident, col.name = "orig.ident")
  sobj <- AddMetaData(sobj, sampleInfo$tissue, col.name = "tissue_type")
  sobj <- AddMetaData(sobj, sampleInfo$woodchuck, col.name = "woodchuck")
  sobj <- AddMetaData(sobj, sampleInfo$sex, col.name = "sex")
  sobj <- AddMetaData(sobj, sampleInfo$date, col.name = "date") 
} else {
  sobj <- AddMetaData(sobj, paste(woodchuck, sample, sep = "_"), col.name = "orig.ident")
  sobj <- AddMetaData(sobj, woodchuck, col.name = "woodchuck")
  sobj <- AddMetaData(sobj, "Female", col.name = "sex")
  sobj <- AddMetaData(sobj, "February", col.name = "date")
}


# Store mitochondrial percentage in object meta data
sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")

# Visualize
VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE)

# # Alternatively, try soupx
# sc <- SoupX::load10X("~/Dropbox/Zoe/scf_version/analysis/multiome/L192_healthy/L192_healthy/outs")
# # For single-nuc
# clusters <- read.csv("~/Dropbox/Zoe/scf_version/analysis/multiome/L192_healthy/L192_healthy/outs/analysis/clustering/gex/graphclust/clusters.csv")
# clusters <- clusters$Cluster
# sc <- SoupX::setClusters(sc, clusters)
# # Can get auto contamination fraction
# sc <- SoupX::autoEstCont(sc, forceAccept = TRUE)
# # OR set it manually
# sc <- SoupX::setContaminationFraction(sc, 0.3)
# sobj <- SoupX::adjustCounts(sc)
# sobj <- CreateSeuratObject(sobj)
# # Also can just estimate soup
# tod <- Read10X("~/Dropbox/Zoe/scf_version/analysis/L202/TLH/L202_TLH_cellranger/outs/raw_feature_bc_matrix/")
# toc <- Read10X("~/Dropbox/Zoe/scf_version/analysis/L202/TLH/L202_TLH_cellranger/outs/filtered_feature_bc_matrix/")
# sc <- SoupX::SoupChannel(tod,toc)
# # Order the soup profile by est (higher est == more in soup)
# soup <- sc$soupProfile[order(-sc$soupProfile$est),]
# soup[1:10,]

# Do this for raw reads instead:
#sobj <- CreateSeuratObject(counts = sobj.data, min.cells = 3,
#                          min.features = 100)

# # Get super high UMI counts
# lowUMIs <- WhichCells(sobj, expression = nCount_RNA < 80000)
# length(lowUMIs)
# # Only keep these cells
# sobj <- subset(sobj, cells = lowUMIs)
# # Now do for UMIs below/above elbow
# lowUMIs <- WhichCells(sobj, expression = nCount_RNA < 766)
# highUMIs <- WhichCells(sobj, expression = nCount_RNA > 765)
# length(highUMIs)
# sobj <- subset(sobj, cells = highUMIs)
# highSobj <- subset(sobj, cells = highUMIs)

# For filtering cells with high mito
# 50% for healthy and infected TLH, 10% for healthy PBMC, 15% for infected PBMC, 2% for single-nuc
if(sample == "PBMC") {
  mitoCells <- WhichCells(sobj, expression = percent.mt < 10)
} else {
  mitoCells <- WhichCells(sobj, expression = percent.mt < 50) 
}
length(mitoCells)
sobj
# Only keep these cells
sobj <- subset(sobj, cells = mitoCells)

# # Alternatively, look just at cells filtered out with high mito
# highMito <- WhichCells(sobj, expression = percent.mt > 50)
# length(highMito)
# # Only keep these cells
# sobj <- subset(sobj, cells = highMito)

# If I want to get rid of the mito genes now
#mitoGenes <- grep("^MT-", sobj@assays$RNA@data@Dimnames[[1]], value = TRUE)
noMitoGenes <- sobj@assays$RNA@data@Dimnames[[1]][!grepl("^MT-", sobj@assays$RNA@data@Dimnames[[1]])]
sobj <- subset(sobj, features = noMitoGenes)

# Do the same for ribosomal genes
sobj <- PercentageFeatureSet(sobj, pattern = "^RPL|^Rpl|^RPS|^Rps", col.name = "percent.ribo")

# Run sctransform
sobj <- SCTransform(sobj) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>% 
  RunUMAP(dims = 1:30) %>%
  RunTSNE()

# Visualization: https://satijalab.org/seurat/articles/visualization_vignette.html

# Perform clustering with scClustViz
scSeurat <- sobj
DE_bw_clust <- TRUE
FDRthresh <- 0.01
seurat_resolution <- 0
sCVdata_list <- list()

################## - SCT assay

while(DE_bw_clust) {
  seurat_resolution <- seurat_resolution + 0.2
  # ^ iteratively incrementing resolution parameter
  
  scSeurat <- FindClusters(scSeurat,
                           resolution=seurat_resolution,
                           print.output=F)
  message(" ")
  message("------------------------------------------------------")
  message(paste0("--------  res.",seurat_resolution," with ",
                 length(levels(Idents(scSeurat)))," clusters --------"))
  message("------------------------------------------------------")
  curr_sCVdata <- CalcSCV(inD=scSeurat,
                          cl=Idents(scSeurat),
                          # ^ your most recent clustering results get stored in the Seurat "ident" slot
                          assayType = "SCT",
                          exponent=exp(1),
                          pseudocount=1,
                          DRthresh=0.1,
                          DRforClust="pca",
                          calcSil=T,
                          calcDEvsRest=T,
                          calcDEcombn=T)
  
  DE_bw_NN <- sapply(DEneighb(curr_sCVdata,FDRthresh),nrow)
  # ^ counts # of DE genes between neighbouring clusters at your selected FDR threshold
  
  if (min(DE_bw_NN) < 1 | seurat_resolution > 2) { DE_bw_clust <- FALSE }
  # ^ If no DE genes between nearest neighbours, don't loop again.
  
  sCVdata_list[[paste0("res.",seurat_resolution)]] <- curr_sCVdata
}

# Now for DropletQC
# Calculate nuclear fractions
nf1 <- nuclear_fraction_tags(
  outs = paste("~/temp_cellranger_outputs/", woodchuck, "/", sample, "/", sep = ""),
  tiles = 1, cores = 1, verbose = FALSE)
clusteringRes <- "SCT_snn_res.0.6"
Idents(scSeurat) <- clusteringRes

# Combine Seurat metadata with nuclear fraction from only the cells in the Seurat object
nuclear_fraction <- nf1[row.names(scSeurat@meta.data),]
fullDF <- cbind(scSeurat@meta.data, nuclear_fraction)
# Get data frame with the nuclear fraction in the first column and umi counts in
# the second; don't need to filter cells by cell type at this point
cells.nf.umi <- data.frame(nf = nuclear_fraction,
                           umi = fullDF$nCount_SCT)
row.names(cells.nf.umi) <- row.names(fullDF)

# Run identify_empty_drops
# Note: damaged cells in bottom right of plot
cells.ed <- identify_empty_drops(nf_umi = cells.nf.umi, include_plot = TRUE)
head(cells.ed)
table(cells.ed$cell_status)

# Identify damaged cells
# Add cell type column to calculate damaged cells
cells.ed$cell_type <- as.character(Idents(scSeurat))
head(cells.ed)
# Identify damaged cells
cells.ed.dc <- identify_damaged_cells(cells.ed, verbose = FALSE, output_plots = TRUE)
table(cells.ed.dc[[1]]$cell_status)
# Which cells are damaged
table(cells.ed.dc[[1]]$cell_type[which(cells.ed.dc[[1]]$cell_status == "damaged_cell")])

# Save DropletQC results as meta data
scSeurat <- AddMetaData(scSeurat, nuclear_fraction, col.name = "nuclear_fraction")
scSeurat <- AddMetaData(scSeurat, cells.ed.dc[[1]]$cell_status, col.name = "dropletQC_result")
#scSeurat <- AddMetaData(scSeurat, cells.ed$cell_status, col.name = "dropletQC_result")

DimPlot(scSeurat, group.by = c("dropletQC_result"), label = TRUE)
FeaturePlot(scSeurat, features = "nuclear_fraction")

boxplot(log10(scSeurat$nCount_RNA) ~ scSeurat$dropletQC_result,
        ylab = "Log10 of nCount_RNA",
        xlab = "DropletQC Result")

# Label with SCINA
# Load GMT
markers <- msigdbi::read.gmt('~/Dropbox/Zoe/scf_version/analysis/pathway_files/woodchuck_cell_markers.gmt')
# Make expression matrix from Seurat object
exprMatrix <- as.matrix(scSeurat@assays$SCT@data)
# Run SCINA - NOTE gene sets need to be 10 or bigger
scina_labels = SCINA::SCINA(exp = exprMatrix, signatures = markers$genesets,
                            rm_overlap = FALSE, allow_unknown = TRUE)
# Add labels as metadata
scSeurat <- Seurat::AddMetaData(scSeurat, scina_labels$cell_labels, col.name = "scina_labels")
# Visualize
DimPlot(scSeurat, group.by = c("scina_labels"), label = TRUE)

# Save scClustViz object
save(sCVdata_list, scSeurat,
     file = paste("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/no_dropletQC/",
                  woodchuck, "/", paste(woodchuck, sample, sep = "_"), "_scClustViz.RData", sep = ""))

# Look at mito distribution for non-empty-droplets
#realCells <- WhichCells(scSeurat, expression = dropletQC_result == "cell")
realCells <- WhichCells(scSeurat, expression = dropletQC_result != "empty_droplet")
realSeurat <- subset(scSeurat, cells = realCells)
VlnPlot(realSeurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, log = TRUE,
        group.by = "orig.ident")
VlnPlot(realSeurat, features = c("percent.ribo", "percent.mt"), ncol = 2,
        group.by = "orig.ident")

# Rerun SCTransform and clustering
realSeurat <- SCTransform(realSeurat) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>% 
  RunUMAP(dims = 1:30) %>%
  RunTSNE()

scSeurat <- realSeurat
DE_bw_clust <- TRUE
FDRthresh <- 0.01
seurat_resolution <- 0
sCVdata_list <- list()

################## - SCT assay

while(DE_bw_clust) {
  seurat_resolution <- seurat_resolution + 0.2
  # ^ iteratively incrementing resolution parameter
  
  scSeurat <- FindClusters(scSeurat,
                           resolution=seurat_resolution,
                           print.output=F)
  message(" ")
  message("------------------------------------------------------")
  message(paste0("--------  res.",seurat_resolution," with ",
                 length(levels(Idents(scSeurat)))," clusters --------"))
  message("------------------------------------------------------")
  curr_sCVdata <- CalcSCV(inD=scSeurat,
                          cl=Idents(scSeurat),
                          # ^ your most recent clustering results get stored in the Seurat "ident" slot
                          assayType = "SCT",
                          exponent=exp(1),
                          pseudocount=1,
                          DRthresh=0.1,
                          DRforClust="pca",
                          calcSil=T,
                          calcDEvsRest=T,
                          calcDEcombn=T)
  
  DE_bw_NN <- sapply(DEneighb(curr_sCVdata,FDRthresh),nrow)
  # ^ counts # of DE genes between neighbouring clusters at your selected FDR threshold
  
  if (min(DE_bw_NN) < 1 | seurat_resolution > 2) { DE_bw_clust <- FALSE }
  # ^ If no DE genes between nearest neighbours, don't loop again.
  
  sCVdata_list[[paste0("res.",seurat_resolution)]] <- curr_sCVdata
}

# Save scClustViz object
save(sCVdata_list, scSeurat,
     file = paste("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/dropletQC_filtered/",
                  woodchuck, "/", sampleInfo$orig.ident, "_scClustViz.RData", sep = ""))

