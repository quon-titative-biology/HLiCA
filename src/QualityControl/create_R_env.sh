### Creating R environment in conda

#conda remove --name R --all
# conda create -n user-base

name='HLiCA_R2'

conda create -n $name -c conda-forge r-base=4.2.3 -y

### Activating R environment
conda activate $name


### installing Mamba for fast downloading of packages in conda
conda install mamba -n base -c conda-forge -y

mamba install -c conda-forge r-essentials --yes

# mamba search r-seurat --channel conda-forge --yes
mamba install r-seurat --channel conda-forge --yes
mamba install r-devtools --channel conda-forge --yes

### Open R and install BiocManager and select a mirror to install the packages from. Use the following

"""
install.packages("remotes")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Update none
BiocManager::install("SingleCellExperiment")
BiocManager::install("GenomeInfoDb")

# Install with conda/mamba
# if (!requireNamespace("devtools", quietly = TRUE))
#     install.packages("devtools")


install.packages("Signac")

BiocManager::install("DropletUtils")
BiocManager::install("scmap")
install.packages('SoupX')
install.packages('ggplot2')
install.packages('optparse')

# For R 4.2.3
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")

# # devtools errors
# # Update CRAN pacakge only
# devtools::install_github("drieslab/Giotto@suite")
# devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
#
"""
