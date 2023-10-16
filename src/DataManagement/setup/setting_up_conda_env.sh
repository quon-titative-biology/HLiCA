"""
"""

# First create and activate conda environment
conda create --name HLiCA4 python=3.9
conda activate HLiCA4

## Update conda
# conda update conda -y
# conda update --all
# conda update -n base -c defaults conda
conda install -c conda-forge mamba -y

mamba install -c conda-forge r-base=4.1.2 -y
# mamba install -c conda-forge r-reticulate">1.0" -y
# mamba install -c conda-forge r-png r-rcurl -y
mamba install -c conda-forge r-seurat -y

# Install seuratdisk
mamba install -c conda-forge libmamba
mamba install -c conda-forge libmambapy

mamba install -c conda-forge r-seuratdisk -y

# mamba install -c r r-stringi
