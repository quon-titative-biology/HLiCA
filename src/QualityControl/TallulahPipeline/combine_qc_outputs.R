library(Seurat)

base='alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc'

# for (dataset in c('DasGupta','Gruen','Henderson','GuilliamsScott','Toronto','Mullen')){
for (dataset in c('DasGupta')){
# for (dataset in c('Henderson')){

  file.names = Sys.glob(file.path(base,dataset,"*",'raw_feature_bc_matrix/qc_output/Cleaned_output_SoupX.rds'))
  # file.names = Sys.glob(file.path('alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc/share_qc',"Toronto","*",'raw_feature_bc_matrix/qc_output/Cleaned_output_SoupX.rds'))
  file.names = sapply(file.names, function(s){dirname(s)})
  # file.names = file.names[c(1,2)]

  i=1

  myseur_list = list()
  empty.sizes = list()

  for (file.name in file.names) {

    print(file.name)

    myseur = readRDS(file.path(file.name,'Cleaned_output_EmptyOnly.rds'))
    rawdata <- Read10X(data.dir = dirname(file.name))
    myseur.soup = readRDS(file.path(file.name,'Cleaned_output_SoupX.rds'))

    # myseur_list[[i]] = myseur.soup

    empty.sizes[[i]] = c(dim(rawdata),dim(myseur),dim(myseur.soup))
    names(empty.sizes[[i]]) = c('Raw # gene', 'Raw # cell', 'EmptyDrop # gene', 'EmptyDrop # cell', 'SoupX # gene','SoupX # cell')
    empty.sizes[[i]][['% Hepatocyte']] = sum(myseur@meta.data$marker_general_labs == "Hepatocyte") / dim(myseur)[[2]]
    empty.sizes[[i]][["EmptyDrop # UMI"]] = sum(myseur@assays$RNA@counts)
    empty.sizes[[i]][["SoupX # UMI"]] = sum(myseur.soup@assays$RNA@counts)
    empty.sizes[[i]][['SoupX % removed']] = 1 - sum(myseur.soup@assays$RNA@counts)/sum(myseur@assays$RNA@counts)

    i=i+1

  }

  # Merge the df stats
  cell.ids = sapply(file.names,function(s){basename(dirname(dirname(s)))})
  df.stat = data.frame(empty.sizes)
  colnames(df.stat) = cell.ids

  csv.file = file.path(base,dataset,'stats.csv')
  write.table(df.stat, file=csv.file,sep=',')

}

# ==============================================================================
#
# ==============================================================================

# Merge and plot
myseur.combined = merge(myseur_list[[1]], myseur_list[-1], add.cell.ids = cell.ids, merge.data = TRUE)



# for (file.name in file.names) {
#
#   print(file.name)
#
#   myseur.soup = readRDS(file.path(file.name,'Cleaned_output_SoupX.rds'))
#
#   myseur_list[[i]] = myseur.soup
#
#   i=i+1
#
# }
#
# # Merge and plot
# cell.ids = sapply(file.names,function(s){basename(dirname(dirname(s)))})
# myseur.combined = merge(myseur_list[[1]], myseur_list[-1], add.cell.ids = cell.ids, merge.data = TRUE)

myseur <- myseur.combined
myseur <- Seurat::ScaleData(myseur);

myseur <- Seurat::FindVariableFeatures(myseur, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(myseur); hvgs <- hvgs[!grepl("^MT-", hvgs)]; # added 22Sept2020
hvgs <- VariableFeatures(myseur); hvgs <- hvgs[!grepl("^Mt-", hvgs)]; # added 22Sept2020
hvgs <- VariableFeatures(myseur); hvgs <- hvgs[!grepl("^mt-", hvgs)]; # added 22Sept2020
VariableFeatures(myseur) <- hvgs;
# PCA
myseur <- Seurat::RunPCA(myseur, features = VariableFeatures(object = myseur))

# Clustering
npcs=20
kNN=20
myseur <- Seurat::FindNeighbors(myseur, dims = 1:npcs)
myseur <- Seurat::FindClusters(myseur, resolution = OPTS$res, k.param=kNN)
# Visualization with TSNE & UMAP
myseur <- Seurat::RunUMAP(myseur, dims = 1:npcs, parallel=FALSE)

##### Make plots ####
agg_coord_by_cluster <- function(coords, clusters) {
  x <- split(seq(nrow(coords)), clusters)
  result <- sapply(x, function(a) apply(coords[a,],2,median))
  return(result)
}

umap_lab_pos <- agg_coord_by_cluster(myseur@reductions$umap@cell.embeddings, myseur@meta.data$seurat_clusters)

file.save = file.path('alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc/share_qc',"Toronto",'umap.png')
png(file.save, width=6, height=6, units="in", res=100)
DimPlot(myseur, reduction="umap", group.by="consistent_labs", pt.size=.1)
dev.off()

orig.ident = myseur@meta.data$cell_ID
myseur@meta.data$orig.ident = sapply(orig.ident,function(s){basename(dirname(dirname(dirname(s))))})

file.save = file.path('alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc/share_qc',"Toronto",'umap_orig_ident.png')
png(file.save, width=12, height=12, units="in", res=100)
# png(file.save)
DimPlot(myseur, reduction="umap", group.by="orig.ident", pt.size=.1)
dev.off()
