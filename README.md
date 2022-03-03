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
Raw sequencing data for this project is available using Qiita: https://qiita.ucsd.edu/ under the study number 10532 (downloading this data will require the creation of an account). After downloading and unzipping the folder, you will find 4 separate FASTQ directories labeled 3823, 3824, 3825, and 3826. Within each folder is an artifact.html file that details which files are forward reads, reverse reads, and barcodes. We need to remove these artifact files for Qiime2's import tool.

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
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path FASTQ/3823 \
  --output-path 3823-single-end-sequences.qza

cd FASTQ/3824
mv Run2_Undetermined_S0_L001_I1_001.fastq.gz barcodes.fastq.gz
mv Run2_Undetermined_S0_L001_R1_001.fastq.gz forward.fastq.gz
mv Run2_Undetermined_S0_L001_R2_001.fastq.gz reverse.fastq.gz
cd ../..
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path FASTQ/3824 \
  --output-path 3824-paired-end-sequences.qza

cd FASTQ/3825
mv Run3_Undetermined_S0_L001_I1_001.fastq.gz barcodes.fastq.gz
mv Run3_Undetermined_S0_L001_R1_001.fastq.gz forward.fastq.gz
mv Run3_Undetermined_S0_L001_R2_001.fastq.gz reverse.fastq.gz
cd ../..
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path FASTQ/3825 \
  --output-path 3825-paired-end-sequences.qza

cd FASTQ/3826
mv Run4_Undetermined_S0_L001_I1_001.fastq.gz barcodes.fastq.gz
mv Run4_Undetermined_S0_L001_R1_001.fastq.gz forward.fastq.gz
mv Run4_Undetermined_S0_L001_R2_001.fastq.gz reverse.fastq.gz
cd ../..
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path FASTQ/3826 \
  --ouput-path 3826-paired-end-sequences.qza
```

Now, we demultiplex the files. Again, procedures for demultiplexing are different depending on if we are working with single-end reads or paired-end reads.
```
qiime demux emp-single \
  --i-seqs 3823-single-end-sequences.qza \
  --m-barcodes-file mapping_files/3823_mapping_file.txt \
  --m-barcodes-column barcode \
  --o-per-sample-sequences 3823_demux.qza \
  --o-error-correction-details 3823_demux-details.qza

qiime demux emp-paired \
  --m-barcodes-file mapping_files/3824_mapping_file.txt \
  --m-barcodes-column barcode \
  --p-rev-comp-mapping-barcodes \
  --i-seqs 3824-paired-end-sequences.qza \
  --o-per-sample-sequences 3824_demux.qza \
  --o-error-correction-details 3824_demux-details.qza

qiime demux emp-paired \
  --m-barcodes-file mapping_files/3825_mapping_file.txt \
  --m-barcodes-column barcode \
  --p-rev-comp-mapping-barcodes \
  --i-seqs 3825-paired-end-sequences.qza \
  --o-per-sample-sequences 3825_demux.qza \
  --o-error-correction-details 3825_demux-details.qza

qiime demux emp-paired \
  --m-barcodes-file mapping_files/3826_mapping_file.txt \
  --m-barcodes-column barcode \
  --p-rev-comp-mapping-barcodes \
  --i-seqs 3826-paired-end-sequences.qza \
  --o-per-sample-sequences 3826_demux.qza \
  --o-error-correction-details 3826_demux-details.qza
```

Now we're going to reorganize all the files.
```
mkdir 3823-preprocessing
mkdir 3824-preprocessing
mkdir 3825-preprocessing
mkdir 3826-preprocessing

mv 3823-single-end-sequences.qza 3823-preprocessing/3823-single-end-sequences.qza
mv 3823_demux-details.qza 3823-preprocessing/3823_demux-details.qza
mv 3823_demux.qza 3823-preprocessing/3823_demux-details.qza

mv 3824-paired-end-sequences.qza 3824-preprocessing/3824-paired-end-sequences.qza
mv 3824_demux-details.qza 3824-preprocessing/3824_demux-details.qza
mv 3824_demux.qza 3824-preprocessing/3824_demux-details.qza

mv 3825-paired-end-sequences.qza 3825-preprocessing/3825-paired-end-sequences.qza
mv 3825_demux-details.qza 3825-preprocessing/3825_demux-details.qza
mv 3825_demux.qza 3825-preprocessing/3825_demux-details.qza

mv 3826-paired-end-sequences.qza 3826-preprocessing/3826-paired-end-sequences.qza
mv 3826_demux-details.qza 3826-preprocessing/3826_demux-details.qza
mv 3826_demux.qza 3826-preprocessing/3826_demux-details.qza
```

Now we do some DADA2 denoising.
```
qiime dada2 denoise-single \
  --i-demultiplexed-seqs 3823-preprocessing/3823_demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --o-representative-sequences 3823-preprocessing/3823-rep-seqs.qza \
  --o-table 3823-preprocessing/3823-table.qza \
  --o-denoising-stats 3823-preprocessing/3823-stats.qza
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 3824-preprocessing/3824_demux.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-representative-sequences 3824-preprocessing/3824-rep-seqs.qza \
  --o-table 3824-preprocessing/3824-table.qza \
  --o-denoising-stats 3824-preprocessing/3824-stats.qza
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 3825-preprocessing/3825_demux.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-representative-sequences 3825-preprocessing/3825-rep-seqs.qza \
  --o-table 3825-preprocessing/3825-table.qza \
  --o-denoising-stats 3825-preprocessing/3825-stats.qza
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 3826-preprocessing/3826_demux.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-representative-sequences 3826-preprocessing/3826-rep-seqs.qza \
  --o-table 3826-preprocessing/3826-table.qza \
  --o-denoising-stats 3826-preprocessing/3826-stats.qza
```

Now that we've completed denoising, we need to merge the tables and the rep-seqs together into a single file. Qiime has built-in commands for the table.qza files and the rep-seq.files. Mapping files need to be handled separately using UNIX commands.
```
qiime feature-table merge \
  --i-tables 3823-preprocessing/3823-table.qza \
  --i-tables 3824-preprocessing/3824-table.qza \
  --i-tables 3825-preprocessing/3825-table.qza \
  --i-tables 3826-preprocessing/3826-table.qza \
  --o-merged-table merged_table.qza
  
qiime feature-table merge-seqs \
  --i-data 3823-preprocessing/3823-rep-seqs.qza \
  --i-data 3824-preprocessing/3824-rep-seqs.qza \
  --i-data 3825-preprocessing/3825-rep-seqs.qza \
  --i-data 3826-preprocessing/3826-rep-seqs \
  --o-merged-data merged_rep-seqs.qza

tail -n +2 mapping_files/3824_mapping_file.txt > mapping_file/cat_3824_mapping_file.txt
tail -n +2 mapping_files/3825_mapping_file.txt > mapping_file/cat_3825_mapping_file.txt
tail -n +2 mapping_files/3826_mapping_file.txt > mapping_file/cat_3826_mapping_file.txt

cat mapping_files/3823_mapping_file.txt \
  mapping_files/cat_3824_mapping_file.txt \
  mapping_files/cat_3825_mapping_file.txt \
  mapping_files/cat_3826_mapping_file.txt \
  > merged_mapping_file.txt
```

Now we can begin our standard qiime analysis.
```
qiime feature-table summarize \
  --i-table merged_table.qza \
  --o-visualization merged_table.qzv \
  --m-sample-metadata-file merged_mapping_file.txt
```

Manual inspection of the feature table reveals a natural sequence depth cutoff at 5892.
![image](https://user-images.githubusercontent.com/57808677/155183751-276b186e-998d-473e-972f-ac8b0bd99512.png)

```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences merged_rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree-qza \
  --i-table merged_table.qza \
  --p-sampling-depth 5892 \
  --m-metadata-file merged_mapping_file.txt \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file merged_mapping_file.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file merged_mapping_file.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv
```

Run the MergeMappingFile.R in my local Capstone folder to modify the metadata text file format.
This will write a new manifest file with extended metadata information. The row "description" contains the patientID followed by a period followed by the time point when the stool sample was collected. These values allow us to identify the alpha diversity of each patient at any given time.
```
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree-qza \
  --i-table merged_table.qza \
  --p-sampling-depth 5892 \
  --m-metadata-file mapping_file_full_metadata.tsv \
  --output-dir core-metrics-results-extended

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-extended/faith_pd_vector.qza \
  --m-metadata-file mapping_file_full_metadata.tsv \
  --o-visualization core-metrics-results-extended/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-extended/evenness_vector.qza \
  --m-metadata-file mapping_file_full_metadata.tsv \
  --o-visualization core-metrics-results-extended/evenness-group-significance.qzv
```

From here I manually unzipped and extracted the "data/alpha-diversity.tsv" file from "core-metrics-results-extended/faith_pd_vector.qza".

I manually curated "clinical_metadata.csv" from supplemental metadata in Kang et al., 2017 [1].

Then I ran the clinical_metadata_analysis.R script to create a linear model between alpha diversity of other metadata metrics. The ggpairs function created clinical_multicollinearity.pdf

## Usage:

## References
[1] Kang, D. W., Adams, J. B., Gregory, A. C., Borody, T., Chittick, L., Fasano, A., Khoruts, A., Geis, E., 
Maldonado, J., McDonough-Means, S., Pollard, E. L., Roux, S., Sadowsky, M. J., Lipson, K. S., Sullivan, M. B., Caporaso, J. G., & Krajmalnik-Brown, R. (2017). Microbiota Transfer Therapy alters gut ecosystem and improves gastrointestinal and autism symptoms: an open-label study. Microbiome, 5(1), 10. https://doi.org/10.1186/s40168-016-0225-7



