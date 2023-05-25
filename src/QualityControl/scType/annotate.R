library('dplyr')
library('Seurat')
library('HGNChelper')
library('openxlsx')

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./alignment/ref_GRCh38p13_gencode_v42/submissions-czi004liv/gruen_2023/MA_2022_2/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
