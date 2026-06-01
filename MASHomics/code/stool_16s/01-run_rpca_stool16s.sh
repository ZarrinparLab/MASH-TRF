#to run RPCA you need to be in environment qiime2-2021.4

#NASH (12wk,NASH) FAFT
path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/
mtb_file=$path/data/stool_16s/mashstool16s_preprocessed_20211020_ID_13785_gg2/NASH_NASH_FAFT_taxonomy_filtered.gg2.asv.counts.qza
md_file=$path/data/stool_16s/s16s_metadata_cln_addmoreNASHcat.txt
output_path=$path/results/stool_16s/s16s_rpca_results_NASH_12wk_FAFT_gg2
comparison=condition
comparison2=collection_timepoint
comparison3=NASH_category
comparison4=fibrosis_stage_new
comparison5=steatosis_grade_new

#pre-NASH (8wk,TRF) FAFT
#path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/
#mtb_file=$path/data/stool_16s/mashstool16s_preprocessed_20211020_ID_13785_gg2/NASH_TRF_FAFT_taxonomy_filtered.gg2.asv.counts.qza
#md_file=$path/data/stool_16s/s16s_metadata_cln_addmoreNASHcat.txt
#output_path=$path/results/stool_16s/s16s_rpca_results_NASH_8wk_FAFT_gg2
#comparison=condition
#comparison2=collection_timepoint
#comparison3=NASH_category
#comparison4=fibrosis_stage_new
#comparison5=steatosis_grade_new

#all NASH
#path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/
#mtb_file=$path/data/stool_16s/mashstool16s_preprocessed_20211020_ID_13785_gg2/NASH_taxonomy_filtered.gg2.asv.counts.qza
#md_file=$path/data/stool_16s/s16s_metadata_cln_addmoreNASHcat.txt
#output_path=$path/results/stool_16s/s16s_rpca_results_NASH_all_gg2
#comparison=condition
#comparison2=collection_timepoint
#comparison3=NASH_category_new
#comparison4=fibrosis_stage_new
#comparison5=steatosis_grade_new

qiime deicode rpca \
  --i-table $mtb_file \
  --p-min-feature-count 0 \
  --p-min-sample-count 0 \
  --o-biplot $output_path/ordination.qza \
  --o-distance-matrix $output_path/distance_matrix.qza

qiime emperor biplot \
  --i-biplot $output_path/ordination.qza \
  --m-sample-metadata-file $md_file \
  --p-ignore-missing-samples \
  --o-visualization $output_path/biplot.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison}-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison2 \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison2}-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison3 \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison3}-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison4 \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison4}-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $output_path/distance_matrix.qza \
  --m-metadata-file $md_file \
  --m-metadata-column $comparison5 \
  --p-method permanova \
  --p-pairwise \
  --o-visualization $output_path/${comparison5}-significance.qzv
