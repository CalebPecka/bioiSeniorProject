# bioiSeniorProject
A Short description about your project

## Installation requirements:

# Step 1. Install a conda environment
An Ubuntu Linux environment for Windows 64 was used for the creation of this project.

```
wget "https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh"

bash Anaconda3-2021.11-Linux-x86_64.sh
```

close and re-open terminal

Test installation of conda
```
conda env list
```

If issues persist, see additional documentation: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

# Step 2. Install qiime2 requirements
wget https://data.qiime2.org/distro/core/qiime2-2021.11-py38-linux-conda.yml

conda env create -n qiime2-2021.11 --file qiime2-2021.11-py38-linux-conda.yml

rm qiime2-2021.11-py38-linux-conda.yml

conda activate qiime2-2021.11

If issues persist, see additional documentation on qiime's website: https://docs.qiime2.org/2021.11/install/native/#install-qiime-2-within-a-conda-environment

# Data Access
Raw sequencing data for this project is available using Qiita: https://qiita.ucsd.edu/ under the study number 10532. After downloading and unzipping the folder (I outputted the contents to a directory called "FASTQ"), you will find 4 separate FASTQ directories labeled 3823, 3824, 3825, and 3826. Within each folder is an artifact.html file that details which files are forward reads, reverse reads, and barcodes. We need to remove these artifact files for Qiime2's import tool.

```
rm -rf FASTQ/3823/artifact_3823.html

rm -rf FASTQ/3824/artifact_3824.html

rm -rf FASTQ/3825/artifact_3825.html

rm -rf FASTQ/3826/artifact_3826.html
```

Qiime2's import tool also requires that we relabel all of the files to a specific convention, which can be followed with the steps below.

```
cd FASTQ/3825

mv Run3_Undetermined_S0_L001_I1_001.fastq.gz barcodes.fastq.gz

mv Run3_Undetermined_S0_L001_R1_001.fastq.gz forward.fastq.gz

mv Run3_Undetermined_S0_L001_R2_001.fastq.gz reverse.fastq.gz

cd ../..

qiime tools import --type EMPPairedEndSequences --input-path FASTQ/3825 --output-path 3825-paired-end-sequences.qza

cd FASTQ/3826

mv Run4_Undetermined_S0_L001_I1_001.fastq.gz barcodes.fastq.gz

mv Run4_Undetermined_S0_L001_R1_001.fastq.gz forward.fastq.gz

mv Run4_Undetermined_S0_L001_R2_001.fastq.gz reverse.fastq.gz

cd ../..

qiime tools import --type EMPPairedEndSequences --input-path FASTQ/3826 --ouput-path 3826-paired-end-sequences.qza
```

## Usage:

## Notes from Kang et al.
"Sequence quality control and demultiplexing using QIIME’s split_libraries_fastq.py with default parameters was performed as described in Bokulich et al. [46] on a per-run basis. The sequences were combined across runs by merging the resulting files using the cat Unix command, and sequences were clustered into operational taxonomic units (OTUs) at sequence similarities of 100 and 97%. "
