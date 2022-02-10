# bioiSeniorProject
A Short description about your project

## Installation requirements:

# Step 1. Install a conda environment
An Ubuntu Linux environment for Windows 64 was used for the creation of this project.

wget "https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh"
bash Anaconda3-2021.11-Linux-x86_64.sh
close and re-open terminal
conda env list (to test)
If issues persist, see additional documentation: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

# Step 2. Install qiime2 requirements
wget https://data.qiime2.org/distro/core/qiime2-2021.11-py38-linux-conda.yml
conda env create -n qiime2-2021.11 --file qiime2-2021.11-py38-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2021.11-py38-linux-conda.yml

conda activate qiime2-2021.11

If issues persist, see additional documentation on qiime's website: https://docs.qiime2.org/2021.11/install/native/#install-qiime-2-within-a-conda-environment

## Usage:

## Notes from Kang et al.
"Sequence quality control and demultiplexing using QIIME’s split_libraries_fastq.py with default parameters was performed as described in Bokulich et al. [46] on a per-run basis. The sequences were combined across runs by merging the resulting files using the cat Unix command, and sequences were clustered into operational taxonomic units (OTUs) at sequence similarities of 100 and 97%. "
