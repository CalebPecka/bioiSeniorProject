# bioiSeniorProject
A Short description about your project

## Installation requirements:

# Step 1. Install a conda environment
An Ubuntu Linux environment for Windows 64 was used for the creation of this project.

```
wget "https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh"
bash Anaconda3-2021.11-Linux-x86_64.sh
```

Close and re-open terminal

Test installation of conda
```
conda env list
```

If issues persist, see additional documentation: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

# Step 2. Install qiime2 requirements
```
wget https://data.qiime2.org/distro/core/qiime2-2021.11-py38-linux-conda.yml
conda env create -n qiime2-2021.11 --file qiime2-2021.11-py38-linux-conda.yml
rm qiime2-2021.11-py38-linux-conda.yml
conda activate qiime2-2021.11
```

If issues persist, see additional documentation on qiime's website: https://docs.qiime2.org/2021.11/install/native/#install-qiime-2-within-a-conda-environment

# Data Access
Raw sequencing data for this project is available using Qiita: https://qiita.ucsd.edu/ under the study number 10532. After downloading and unzipping the folder, you will find 4 separate FASTQ directories labeled 3823, 3824, 3825, and 3826. Within each folder is an artifact.html file that details which files are forward reads, reverse reads, and barcodes. We need to remove these artifact files for Qiime2's import tool.

```
rm -rf FASTQ/3823/artifact_3823.html
rm -rf FASTQ/3824/artifact_3824.html
rm -rf FASTQ/3825/artifact_3825.html
rm -rf FASTQ/3826/artifact_3826.html
```

Qiime2's import tool also requires that we relabel all of the files to a specific convention, which can be followed with the steps below. It is important to note that FASTQ reads for directories 3824, 3825, and 3826 use paired end sequence reads, while directory 3823 uses single end sequence reads. Qiime2 commands up until the demultiplexing step must be carried out separately.

**Notes from Kang et al. 2017 [1]**
"Sequence quality control and demultiplexing using QIIME’s split_libraries_fastq.py with default parameters was performed as described in Bokulich et al. [46] on a per-run basis. The sequences were combined across runs by merging the resulting files using the cat Unix command, and sequences were clustered into operational taxonomic units (OTUs) at sequence similarities of 100 and 97%."

```
cd FASTQ/3823
mv Undetermined_S0_L001_I1_001.fastq.gz barcodes.fastq.gz
mv Undetermined_S0_L001_R1_001.fastq.gz sequences.fastq.gz
cd ../..
qiime tools import --type EMPSingleEndSequences --input-path FASTQ/3823 --output-path 3823-single-end-sequences.qza

cd FASTQ/3824
mv Run2_Undetermined_S0_L001_I1_001.fastq.gz barcodes.fastq.gz
mv Run2_Undetermined_S0_L001_R1_001.fastq.gz forward.fastq.gz
mv Run2_Undetermined_S0_L001_R2_001.fastq.gz reverse.fastq.gz
cd ../..
qiime tools import --type EMPPairedEndSequences --input-path FASTQ/3824 --output-path 3824-paired-end-sequences.qza

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

Now, we demultiplex the files. Again, procedures for demultiplexing are different depending on if we are working with single-end reads or paired-end reads.
```
qiime demux emp-single --i-seqs 3823-single-end-sequences.qza --m-barcodes-file mmaping_files/3823_mapping_file.txt --m-barcodes-column barcode --o-per-sample-sequences 3823_demux.qza --o-error-correction-details 3823_demux-details.qza
```
## Usage:

## References
[1] Kang, D. W., Adams, J. B., Gregory, A. C., Borody, T., Chittick, L., Fasano, A., Khoruts, A., Geis, E., 
Maldonado, J., McDonough-Means, S., Pollard, E. L., Roux, S., Sadowsky, M. J., Lipson, K. S., Sullivan, M. B., Caporaso, J. G., & Krajmalnik-Brown, R. (2017). Microbiota Transfer Therapy alters gut ecosystem and improves gastrointestinal and autism symptoms: an open-label study. Microbiome, 5(1), 10. https://doi.org/10.1186/s40168-016-0225-7



