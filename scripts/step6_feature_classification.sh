cd ..

qiime tools import --type 'FeatureData[Sequence]' \
  --input-path bac120_ssu_reps_r202.fna \
  --output-path bac120_ssu_reps_r202.qza
  
qiime tools import --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path bac120_taxonomy_r202.tsv \
  --output-path bac120_taxonomy_r202.qza

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