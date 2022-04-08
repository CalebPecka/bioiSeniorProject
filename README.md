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
  --type MultiplexedSingleEndBarcodeInSequence \
  --input-path FASTQ/3823/sequences.fastq.gz \
  --output-path 3823_multiplexed_sequences.qza

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
qiime cutadapt demux-single \
  --i-seqs 3823_multiplexed_seqs.qza \
  --m-barcodes-file mapping_files/3823_mapping_file.txt \
  --m-barcodes-column barcode \
  --p-error-rate 0 \
  --o-per-sample-sequences 3823_demux.qza \
  --o-untrimmed-sequences 3823_untrimmed.qza

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

mv 3823_untrimmed.qza 3823-preprocessing/3823_untrimmed.qza
mv 3823_multiplexed_seqs.qza 3823-preprocessing/3823_multiplexed_seqs.qza
mv 3823_demux.qza 3823-preprocessing/3823_demux.qza

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
  --i-data 3826-preprocessing/3826-rep-seqs.qza \
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

Manual inspection of the feature table reveals that many samples have poor feature coverage. We want to preserve as many samples as reasonably possible since we need week18 results AND week0 results to make meaningful comparisons. For this reason, a low sampling depth of 1525 was chosen for this experiment.
![SequencingDepth](https://user-images.githubusercontent.com/57808677/162447178-6bf9a34c-12ab-422d-bd45-0fb5af44d792.PNG)

```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences merged_rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

Run the MergeMappingFile.R in my local Capstone folder to modify the metadata text file format.
This will write a new manifest file with extended metadata information. The row "description" contains the patientID followed by a period followed by the time point when the stool sample was collected. These values allow us to identify the alpha diversity of each patient at any given time.
```
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table merged_table.qza \
  --p-sampling-depth 1525 \
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
```
     faith_pd --- week           : -0.138*
          bmi --- height         :  0.637***
          bmi --- weight         :  0.905***
          bmi --- age            :  0.646***
          bmi --- weight.pounds  :  0.904***
          bmi --- PGI.improve    :  0.283***
       height --- weight         :  0.890***
       height --- age            :  0.877***
       height --- weight.pounds  :  0.890***
       height --- ABC.change     : -0.183*
       height --- ABC.change.1   :  0.223*
       height --- SRS.change     :  0.190*
       weight --- age            :  0.810***
       weight --- weight.pounds  :  1.000***
       weight --- PGI.improve    :  0.154*
          age --- weight.pounds  :  0.800***
          age --- ABC.change     : -0.244**
weight.pounds --- PGI.improve    :  0.155*
```
 
There were major correlations *** between all metrics of improvement (ABC versus SRS versus PGI R).
There were no correlations between any of these metrics and alpha diversity. Week is the only statistically significant correlation with alpha diversity.
  
Based on these results: bmi, height, and weight.pounds were removed.

![FaithPD_LinearModel](https://user-images.githubusercontent.com/57808677/162479785-f838a2a8-e178-4ec8-876e-f6b6b2f056e2.PNG)

Age, weight, and gender are all having a significant impact on the microbiome alpha diversity.

Given there is no multicollinearity between patient faith_pd and clinical symptoms, faith_pd is a poor measurement of microbiome recovery.

Violin plots were used to show that, on average, faith_pd alpha diversity is decreasing after FMT, becoming more similar to neurotypical individuals. Wilcox unpaired tests confirmed this observation.

![FaithPD_violin_week0ASD](https://user-images.githubusercontent.com/57808677/162483601-aeedd775-73f4-49a8-9d60-711c294b4a66.png)

![FaithPD_violin_week18ASD](https://user-images.githubusercontent.com/57808677/162483612-4b2725cb-f035-43b5-96c3-6cd4d8775833.png)

More linear models at the end of clinical_metadata_analysis.R were produced, separating neurotypical individuals and autism individuals based on week 0 versus week 18. In general, these models were noisy, and the F statistics were not significant. These models were based on age, weight, and gender. Individuals vaues were statistically significant in week 18.

All of these observations were made using faith's phylogenetic diversity score. The literature I've read that claims alpha diversity increases is based on shannon's diversity. "Shannon alpha diversity is sensitive to both the richness (total number of species in the community) and the evenness (relative abundance of different species). Faith's phylogenetic diversity represents the number of phylogenetic tree-units within a sample." (https://www.researchgate.net/figure/Boxplots-showing-distribution-of-Shannon-and-Faiths-phylogenetic-alpha-diversity_fig2_315813036#:~:text=Shannon%20alpha%20diversity%20is%20sensitive,tree%2Dunits%20within%20a%20sample.)

The next logical step is to reperform this analysis using shannon alpha diversity indices. To do this, I ran the Rscript for "shannon_alpha_diversity_analysis.R".

"alpha-diversity.tsv" was manually extracted from "core-metrics-results-extended/shannon_vector.qza" and renamed to "alpha-diversity-shannon.tsv".

Multicollinearity analysis from "graphics/clinical_multicollinearity_shannon.pdf" shows similar results as faith_pd, **except all clinical recovery outcomes are correlated with shannon alpha diversity.** I also note that in both faith_pd and shannon diversity, **measurements of clinical improvement were sometimes significantly correlated with age, weight, and/or height.**

Based on violin plots, shannon entropy is increasing in autism patients from week 0 to week 18, stabilizing **higher** than the expected distribution in neurotypical individuals. This is especially prominent if we submit neurotypical samples to only include week 0. This result is supportive of the claim that FMT is not providing the patients with a more evenly distributed microbiome, it is always introducing higher diversity of bacteria.

![Shannon_violin_week0ASD](https://user-images.githubusercontent.com/57808677/162488132-ee34e094-974b-41a3-ab1d-04f0da13a3d7.png)

![Shannon_violin_week18ASD](https://user-images.githubusercontent.com/57808677/162488139-472e316c-4e11-456a-b9c8-c5dd879ceace.png)

![Shannon_violin_week0ASD_week0neurotypical](https://user-images.githubusercontent.com/57808677/162488158-16b7bf63-fb08-4f22-a630-e428b8b5cb17.png)

![Shannon_violin_week18ASD_week0neurotypical](https://user-images.githubusercontent.com/57808677/162488168-16d0930b-676a-4ddf-96d9-11827f10949a.png)




Linear model of week 0 neurotypical individuals:

![LM_neurotypical_week0](https://user-images.githubusercontent.com/57808677/162489421-cfe8d8c8-5231-41a7-b8f8-46c51d435ef7.png)

```
Residuals:
    Min      1Q  Median      3Q     Max 
-1.7325 -0.5592  0.1107  0.7350  1.4542 

Coefficients:
                Estimate Std. Error t value Pr(>|t|)   
(Intercept)      4.44030    1.33956   3.315  0.00346 **
age             -0.20574    0.16321  -1.261  0.22195   
genderM          0.61085    0.80402   0.760  0.45627   
weight..pounds.  0.02098    0.01343   1.562  0.13393   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.064 on 20 degrees of freedom
Multiple R-squared:  0.1547,	Adjusted R-squared:  0.02789 
F-statistic:  1.22 on 3 and 20 DF,  p-value: 0.3283
```

Linear model of week 18 neurotypical individuals:

![LM_neurotypical_week18](https://user-images.githubusercontent.com/57808677/162489447-729279a3-05fa-48db-8b85-729381012cf5.png)

```
Residuals:
     Min       1Q   Median       3Q      Max 
-2.37082 -0.65415  0.00453  0.87375  2.10216 

Coefficients:
                Estimate Std. Error t value Pr(>|t|)   
(Intercept)      4.42435    1.28519   3.443  0.00212 **
age             -0.21672    0.15701  -1.380  0.18022   
genderM          0.99254    0.64518   1.538  0.13704   
weight..pounds.  0.02039    0.01592   1.281  0.21250   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.176 on 24 degrees of freedom
Multiple R-squared:  0.181,	Adjusted R-squared:  0.07866 
F-statistic: 1.768 on 3 and 24 DF,  p-value: 0.1801
```

Linear model of week 0 ASD individuals:

![LM_ASD_week0](https://user-images.githubusercontent.com/57808677/162489502-1c60a7e7-2577-4305-8e98-701075464fc2.png)

```
Residuals:
     Min       1Q   Median       3Q      Max 
-0.84595 -0.32739 -0.05381  0.23401  0.86596 

Coefficients:
                 Estimate Std. Error t value Pr(>|t|)   
(Intercept)      3.242969   0.845876   3.834  0.00499 **
age              0.297086   0.133748   2.221  0.05707 . 
weight..pounds. -0.017166   0.009508  -1.806  0.10863   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.5667 on 8 degrees of freedom
Multiple R-squared:  0.3854,	Adjusted R-squared:  0.2318 
F-statistic: 2.508 on 2 and 8 DF,  p-value: 0.1427
```

Linear model of week 18 ASD individuals:

![LM_ASD_week18](https://user-images.githubusercontent.com/57808677/162489529-5a7a9476-e30d-4e80-8927-a56f777fb006.png)

```
Residuals:
    Min      1Q  Median      3Q     Max 
-2.3974 -0.1181  0.1110  0.4275  1.2149 

Coefficients:
                 Estimate Std. Error t value Pr(>|t|)    
(Intercept)      4.205184   0.893379   4.707 0.000135 ***
age              0.130116   0.123581   1.053 0.304951    
weight..pounds. -0.005196   0.008158  -0.637 0.531358    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.8598 on 20 degrees of freedom
Multiple R-squared:  0.05812,	Adjusted R-squared:  -0.03607 
F-statistic: 0.6171 on 2 and 20 DF,  p-value: 0.5495
```












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



