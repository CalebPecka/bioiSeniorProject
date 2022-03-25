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

tail -n +2 mapping_files/3824_mapping_file.txt > mapping_files/cat_3824_mapping_file.txt
tail -n +2 mapping_files/3825_mapping_file.txt > mapping_files/cat_3825_mapping_file.txt
tail -n +2 mapping_files/3826_mapping_file.txt > mapping_files/cat_3826_mapping_file.txt

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

Significant correlations include:

  week AND faith_pd -> -0.123*
  
  height AND bmi -> 0.653***
  
  weight AND bmi -> 0.915***
  
  age AND bmi -> 0.686***
  
  weight.pounds AND bmi -> 0.915***
  
  end of treatment PGI R improvement AND bmi -> 0.226**
  
  weight AND height -> 0.888***
  
  age AND height -> 0.890***
  
  weight.pounds AND height -> 0.888***
  
  ABC change AND height -> 0.234**
  
  SRS change AND height -> 0.232**
  
  age AND weight -> 0.831***
  
  weight.pounds AND weight -> 1.000***
  
  weight.pounds AND age -> 0.810***
  
  ABC change AND age -> -0.209**
  
  
  major correlations *** between all metrics of improvement (ABC versus SRS versus PGI R)
  There were no correlations between any of these metrics and alpha diversity
  
Based on these results: bmi, height, and weight.pounds were removed.
One linear model was created without metrics of improvment (ABC versus SRC versus PGI R), and another model included those metrics.

Neither model could explain variation in alpha diversity, signifying that alpha diversity may not be a useful metric for measuring patient recovery.
Here is a graphic of the linear model of age_versus_alpha_diversity in neurotypical individuals
![image](https://github.com/CalebPecka/bioiSeniorProject/blob/main/graphics/neurotypical_age_versus_alpha_diversity.pdf)

And here is the same graphic of the linear model in autism individuals
![image](https://github.com/CalebPecka/bioiSeniorProject/blob/main/graphics/autism_age_versus_alpha_diversity.pdf)

More linear models at the end of clinical_metadata_analysis.R were produced, separating neurotypical individuals and autism individuals based on week 0 versus week 18. Alpha diversity overall decreased in BOTH experimental groups after the 18 week period. In addition, linear models were created for week0 neurotypical patients and week18 neurotypical patients, including metadata values for age, weight, and gender. Week0 neurotypical patients did NOT have a significant change in alpha diversity based on age, but week18 patients DID. This same observation is observed in autism patients. Overall, the treatment appears to be lowering alpha diversity, and the reduce in alpha diversity scales relative to age. (Violin plots are in graphics folder as week0 and week18_autism_violin.pdf)

All of these observations were made using faith's phylogenetic diversity score. The literature I've read that claims alpha diversity increases is based on shannon's diversity. "Shannon alpha diversity is sensitive to both the richness (total number of species in the community) and the evenness (relative abundance of different species). Faith's phylogenetic diversity represents the number of phylogenetic tree-units within a sample." (https://www.researchgate.net/figure/Boxplots-showing-distribution-of-Shannon-and-Faiths-phylogenetic-alpha-diversity_fig2_315813036#:~:text=Shannon%20alpha%20diversity%20is%20sensitive,tree%2Dunits%20within%20a%20sample.)

The next logical step is to reperform this analysis using shannon alpha diversity indices. To do this, I ran the Rscript for "shannon_alpha_diversity_analysis.R".

"alpha-diversity.tsv" was manually extracted from "core-metrics-results-extended/shannon_vector.qza" and renamed to "alpha-diversity-shannon.tsv".

Multicollinearity analysis from "graphics/clinical_multicollinearity_shannon.pdf" shows similar results as faith_pd, **except all clinical recovery outcomes are correlated with shannon alpha diversity.**

Based on violin plots, shannon entropy is increasing in autism patients from week 0 to week 18, stabilizing closer to the expected distribution in neurotypical individuals.





Linear model of week0 neurotypical individuals:
Coefficients:
                Estimate Std. Error t value Pr(>|t|)   
(Intercept)      4.89984    1.50509   3.256  0.00466 **
age             -0.29831    0.19603  -1.522  0.14646   
genderM          0.45694    0.86489   0.528  0.60410   
weight..pounds.  0.02867    0.01635   1.754  0.09744 . 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.115 on 17 degrees of freedom
Multiple R-squared:  0.1992,	Adjusted R-squared:  0.05793 
F-statistic:  1.41 on 3 and 17 DF,  p-value: 0.2743




Linear model of week18 neurotypical individuals:
Coefficients:
                Estimate Std. Error t value Pr(>|t|)   
(Intercept)      4.44091    1.29573   3.427   0.0022 **
age             -0.22119    0.15830  -1.397   0.1751   
genderM          1.03042    0.65047   1.584   0.1263   
weight..pounds.  0.02095    0.01605   1.305   0.2042   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.185 on 24 degrees of freedom
Multiple R-squared:  0.1877,	Adjusted R-squared:  0.0862 
F-statistic: 1.849 on 3 and 24 DF,  p-value: 0.1653






Linear model of week0 autism individuals:
Coefficients:
                                      Estimate Std. Error t value Pr(>|t|)    
(Intercept)                           4.692563   0.545385   8.604 0.000136 ***
age                                   0.224746   0.076111   2.953 0.025518 *  
weight..pounds.                      -0.018778   0.003925  -4.784 0.003051 ** 
end.of.treatment.PGI.R..improvement. -0.309521   0.114851  -2.695 0.035815 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.319 on 6 degrees of freedom
Multiple R-squared:  0.8532,	Adjusted R-squared:  0.7798 
F-statistic: 11.62 on 3 and 6 DF,  p-value: 0.00653

Linear regression for the above plot based on age is in "graphics/shannon_diversity_decreases_with_age_in_autism_week0.pdf"
[shannon_diversity_decreases_with_age_in_autism_week0.pdf](https://github.com/CalebPecka/bioiSeniorProject/files/8187427/shannon_diversity_decreases_with_age_in_autism_week0.pdf)






Linear model of week18 autism individuals:
Coefficients:
                                     Estimate Std. Error t value Pr(>|t|)    
(Intercept)                          4.545769   0.878602   5.174 7.63e-05 ***
age                                  0.022494   0.127143   0.177    0.862    
weight..pounds.                      0.001599   0.008572   0.186    0.854    
end.of.treatment.PGI.R..improvement. 0.192399   0.230900   0.833    0.416    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.7671 on 17 degrees of freedom
Multiple R-squared:  0.08139,	Adjusted R-squared:  -0.08072 
F-statistic: 0.5021 on 3 and 17 DF,  p-value: 0.6859




"distance-matrix.tsv" was manually extracted from "core-metrics-results-extended/unweighted_unifrac_distance_matrix.qza".



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

Now we import these files for classification.
```
qiime tools import --type 'FeatureData[Sequence]' \
  --input-path bac120_ssu_reps_r202.fna \
  --output-path bac120_ssu_reps_r202.qza
  
qiime tools import --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path bac120_taxonomy_r202.tsv \
  --output-path bac120_taxonomy_r202.qza
```

Now we run feature classification to determine the taxonomy of our ASV sequences.
```
qiime feature-classifier classify-consensus-blast \
  --i-query merged_rep-seqs.qza \
  --i-reference-reads bac120_ssu_reps_r202.qza \
  --i-reference-taxonomy bac120_taxonomy_r202.qza \
  --p-perc-identity 0.8 --p-min-consensus 0.51 \
  --p-maxaccepts 10 \
  --o-classification taxonomy.qza
  
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table merged_table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file merged_mapping_file.txt \
  --o-visualization taxa-bar-plots.qzv
```

"level-7.csv" was manually extracted from the zipped archive for "taxa-bar-plots.qzv".

Now we're interested in distinguishing individuals whose microbiome recovered versus those who did not recover. 

"distance-matrix.tsv" was manually extracted from the zipped archive for "core-metrics-results-extended/jaccard_distance_matrix.qza".

## Usage:

## References
[1] Kang, D. W., Adams, J. B., Gregory, A. C., Borody, T., Chittick, L., Fasano, A., Khoruts, A., Geis, E., 
Maldonado, J., McDonough-Means, S., Pollard, E. L., Roux, S., Sadowsky, M. J., Lipson, K. S., Sullivan, M. B., Caporaso, J. G., & Krajmalnik-Brown, R. (2017). Microbiota Transfer Therapy alters gut ecosystem and improves gastrointestinal and autism symptoms: an open-label study. Microbiome, 5(1), 10. https://doi.org/10.1186/s40168-016-0225-7



