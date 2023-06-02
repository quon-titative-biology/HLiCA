"""
"""

# First create and activate conda environment
conda create --name hlica
conda activate hlica

## Update conda
# conda update conda -y
# conda update --all
# conda update -n base -c defaults conda
conda install -c conda-forge mamba -y


mamba install -c conda-forge r-base=4.1.2
conda install r-png r-rcurl -y

mamba install -c r r-stringi
