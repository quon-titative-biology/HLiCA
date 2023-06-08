### installing Mamba for fast downloading of packages in conda
`conda install mamba -n base -c conda-forge -y`


### Creating R environment in conda

```
#conda remove --name R --all
conda create -n user-base
mamba create -n R -c conda-forge r-base=4.2.2 -y
```

### Activating R environment
```
conda activate R
mamba install -c conda-forge r-essentials
```

### Open R and install BiocManager and select a mirror to install the packages from. Use the following
```
install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
```
```
install.packages('Seurat')
BiocManager::install("SingleCellExperiment")
remotes::install_github('satijalab/seurat-wrappers')
BiocManager::install("GenomeInfoDb")
install.packages("Signac")
devtools::install_github("drieslab/Giotto@suite")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
```
