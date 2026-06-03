#to run RPCA you need to be in environment qiime2-2021.4

#NASH (12wk,NASH) FAFT
path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/
mtb_file=$path/data/stool_mtb/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk12-FAFT.qza
md_file=$path/data/stool_mtb/NASH_stool_metabolomics_metadata_new.txt
output_path=$path/results/stool_mtb/smtb_rpca_results_NASH_Wk12_FAFT
comparison=ATTRIBUTE_condition
comparison2=ATTRIBUTE_timepoint
comparison3=NASH_category
comparison4=fibrosis_stage_new
comparison5=steatosis_grade_new

#pre-NASH (8wk,TRF) FAFT
#path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/
#mtb_file=$path/data/stool_mtb/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk8-FAFT.qza
#md_file=$path/data/stool_mtb/NASH_stool_metabolomics_metadata_new.txt
#output_path=$path/results/stool_mtb/smtb_rpca_results_NASH_Wk8_FAFT
#comparison=ATTRIBUTE_condition
#comparison2=ATTRIBUTE_timepoint
#comparison3=NASH_category
#comparison4=fibrosis_stage_new
#comparison5=steatosis_grade_new

#all NASH
#path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/MASH-TRF/MASHomics/
#mtb_file=$path/data/stool_mtb/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-rmblnk.qza
#md_file=$path/data/stool_mtb/NASH_stool_metabolomics_metadata_new.txt
#output_path=$path/results/stool_mtb/smtb_rpca_results_NASH_all
#comparison=ATTRIBUTE_condition
#comparison2=ATTRIBUTE_timepoint
#comparison3=NASH_category
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
