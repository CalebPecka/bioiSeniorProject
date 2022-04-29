cd ..

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table merged_table.qza \
  --p-sampling-depth 1525 \
  --m-metadata-file data/mapping_file_full_metadata.tsv \
  --output-dir core-metrics-results-extended

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-extended/faith_pd_vector.qza \
  --m-metadata-file data/mapping_file_full_metadata.tsv \
  --o-visualization core-metrics-results-extended/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-extended/evenness_vector.qza \
  --m-metadata-file data/mapping_file_full_metadata.tsv \
  --o-visualization core-metrics-results-extended/evenness-group-significance.qzv