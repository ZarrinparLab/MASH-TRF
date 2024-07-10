#to run RPCA you need to be in environment qiime2-2021.4

#NASH FT
path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
mtb_file=$path/files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/cecum/data/raw-data-pre-10956-artifact-119356/cecum_NASH_FT_cleaned_119356_reference-hit.qza
md_file=$path/16s/cecum/metadata_cln_addmoreNASHcat.txt
output_path=$path/16s/cecum/rpca_results/rpca_results_NASH_FT
comparison=collection_timepoint
comparison2=NASH_category
comparison3=fibrosis_stage_new
comparison4=steatosis_grade_new

#NASH FA
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/cecum/data/raw-data-pre-10956-artifact-119356/cecum_NASH_FA_cleaned_119356_reference-hit.qza
#md_file=$path/16s/cecum/metadata_cln_addmoreNASHcat.txt
#output_path=$path/16s/cecum/rpca_results/rpca_results_NASH_FA
#comparison=collection_timepoint
#comparison2=NASH_category
#comparison3=fibrosis_stage_new
#comparison4=steatosis_grade_new

#NASH NA
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/cecum/data/raw-data-pre-10956-artifact-119356/cecum_NASH_NA_cleaned_119356_reference-hit.qza
#md_file=$path/16s/cecum/metadata_cln.txt
#output_path=$path/16s/cecum/rpca_results/rpca_results_NASH_NA
#comparison=collection_timepoint

#NASH ZT13
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/cecum/data/raw-data-pre-10956-artifact-119356/cecum_NASH_ZT13_cleaned_119356_reference-hit.qza
#md_file=$path/16s/cecum/metadata_cln_addmoreNASHcat.txt
#output_path=$path/16s/cecum/rpca_results/rpca_results_NASH_ZT13
#comparison=condition
#comparison2=NASH_category
#comparison3=fibrosis_stage_new
#comparison4=steatosis_grade_new

#NASH ZT1
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/cecum/data/raw-data-pre-10956-artifact-119356/cecum_NASH_ZT1_cleaned_119356_reference-hit.qza
#md_file=$path/16s/cecum/metadata_cln_addmoreNASHcat.txt
#output_path=$path/16s/cecum/rpca_results/rpca_results_NASH_ZT1
#comparison=condition
#comparison2=NASH_category
#comparison3=fibrosis_stage_new
#comparison4=steatosis_grade_new

#all NASH
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/cecum/data/raw-data-pre-10956-artifact-119356/cecum_NASH_cleaned_119356_reference-hit.qza
#md_file=$path/16s/cecum/metadata_cln_addmoreNASHcat.txt
#output_path=$path/16s/cecum/rpca_results/rpca_results_NASH_all
#comparison=condition
#comparison2=collection_timepoint
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
