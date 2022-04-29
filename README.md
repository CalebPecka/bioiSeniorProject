# Investigating stagnant clinical outcomes after fecal microbiome transplant in autism spectrum disorder
Research evidence has shown that there is a correlation between the microbial composition of the gut microbiota and function in the brain, a phenomenon called the “gut-brain axis.” The bi-directional correlation between the gut and brain has appeared in unexpected domains and diseases, including gastrointestinal infections that frequently coexist with autism spectrum disorder (ASD) (Ding et al., 2017; Vuong et al., 2017). The function of the gut-brain axis is attributed to bacterial metabolites. Most prominently, metabolites can enter the circulatory system and travel to the brain, or they can be directly transferred through the vagus nerve (Garcia et al., 2021).

ASD is known to have a strong genetic linkage, including a reduced susceptibility to genetic mutations (Iossifov et al., 2015). Studies have also linked the functionality of specific genes involved in ASD to the microbiome composition of the same patients (Sgritta et al., 2019). In an effort to treat ASD symptoms, numerous studies have performed fecal microbiome transplants (FMT) in ASD patients and observed the effects of their recovery. Researchers often focus on the success of these trials, and poorly investigate patients that did not recover after FMT. The goal of this study is to determine several criteria for FMT recovery and investigate possible biological explanations for stagnant ASD recovery outcomes. We hypothesize that variations in FMT donor profiles as well as patient diet may explain variability in critical metabolites associated with ASD FMT recovery.

## Installation requirements:

# Step 1. Install a conda environment

An Ubuntu Linux environment for Windows 64 was used for the creation of this project.

```
wget "https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh"
bash Anaconda3-2021.11-Linux-x86_64.sh
```

Close and re-open your terminal.

Test your installation of conda.

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

# Step 3. Data replication and reproducibility

# Data Access
Raw sequencing data for this project is available using Qiita: https://qiita.ucsd.edu/ under the study number 10532 (downloading this data will require the creation of an account). After downloading and unzipping the folder, you will find 4 separate FASTQ directories labeled 3823, 3824, 3825, and 3826. Within each folder is an artifact.html file that details which files are forward reads, reverse reads, and barcodes. We need to remove these artifact files for Qiime2's import tool.

```
rm -rf FASTQ/3823/artifact_3823.html
rm -rf FASTQ/3824/artifact_3824.html
rm -rf FASTQ/3825/artifact_3825.html
rm -rf FASTQ/3826/artifact_3826.html
```

Qiime2's import tool also requires that we relabel all of the files to a specific convention, which can be followed with the steps below. It is important to note that FASTQ reads for directories 3824, 3825, and 3826 use paired end sequence reads, while directory 3823 uses single end sequence reads. Qiime2 commands up until the demultiplexing step must be carried out separately.

GTDB release202/202.0 was downloaded using the following command.
```
curl -LJO "https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/bac120_taxonomy_r202.tsv.gz"
```

GTDB reference sequences were downloaded using the following command.
```
curl -LJO "https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/genomic_files_reps/bac120_ssu_reps_r202.tar.gz"
```

Then we need to unzip their contents.
```
gunzip bac120_taxonomy_r202.tsv.gz
tar -zxvf bac120_ssu_reps_r202.tar.gz
```

# Usage
The following steps show the order for the analysis of this project:

- Run 'scripts/step1_qiime_preprocessing.sh'
- Run 'scripts/step2_merge_mapping_file.R'
- Run 'scripts/step3_qiime_diversity.sh'
- Manually unzip and extract 'data/alpha-diversity.tsv' from 'core-metrics-results-extended/faith_pd_vector.qza'
- Run 'scripts/step4_faithpd_alpha_diversity_analysis.R'
- Run 'scripts/step5_shannon_alpha_diversity_analysis.R'
- Manually unzip and extract 'data/alpha-diversity.tsv' from 'core-metrics-results-extended/shannon_vector.qza' and rename the file to 'data/alpha-diversity-shannon.tsv'
- Manually unzip and extract 'data/distance-matrix.tsv' from 'core-metrics-results-extended/unweighted_unifrac_distance_matrix.qza'
- Run 'scripts/step6_feature_classification.sh'
- Manually unzip and extract 'level-7.csv' from 'taxa-bar-plots.qzv'
- Manually unzip and extract 'distance-matrix.tsv' from 'core-metrics-results-extended/jaccard_distance_matrix.qza'
- Run 'scripts/step7_differential_abundance_analysis.R'

Comments and plots can be overwritten and saved from 'scripts/step7' to produce different results in different subgroups.

## References
[1] Kang, D. W., Adams, J. B., Gregory, A. C., Borody, T., Chittick, L., Fasano, A., Khoruts, A., Geis, E., 
Maldonado, J., McDonough-Means, S., Pollard, E. L., Roux, S., Sadowsky, M. J., Lipson, K. S., Sullivan, M. B., Caporaso, J. G., & Krajmalnik-Brown, R. (2017). Microbiota Transfer Therapy alters gut ecosystem and improves gastrointestinal and autism symptoms: an open-label study. Microbiome, 5(1), 10. https://doi.org/10.1186/s40168-016-0225-7



