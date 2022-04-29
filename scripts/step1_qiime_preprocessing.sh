cd ..
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

qiime feature-table summarize \
  --i-table merged_table.qza \
  --o-visualization merged_table.qzv \
  --m-sample-metadata-file merged_mapping_file.txt

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences merged_rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza