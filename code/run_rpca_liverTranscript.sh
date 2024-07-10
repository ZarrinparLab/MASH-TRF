#to run RPCA you need to be in environment qiime2-2021.4

#NASH FT
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/transcriptomics/data_files/liver_FT_gene_count_matrix.tsv
#md_file=$path/transcriptomics/liver_metadata_wnashscore_fibste.tsv
#output_path=$path/transcriptomics/rpca_results_liver/rpca_results_NASH_FT
#comparison=timepoint
#comparison2=NASH_category
#comparison3=fibrosis_stage_new
#comparison4=steatosis_grade_new

#NASH FA
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/transcriptomics/data_files/liver_FA_gene_count_matrix.tsv
#md_file=$path/transcriptomics/liver_metadata_wnashscore_fibste.tsv
#output_path=$path/transcriptomics/rpca_results_liver/rpca_results_NASH_FA
#comparison=timepoint
#comparison2=NASH_category
#comparison3=fibrosis_stage_new
#comparison4=steatosis_grade_new

#NASH NA
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/transcriptomics/data_files/liver_NA_gene_count_matrix.tsv
#md_file=$path/transcriptomics/liver_metadata_wnashscore_fibste.tsv
#output_path=$path/transcriptomics/rpca_results_liver/rpca_results_NASH_NA
#comparison=timepoint

#NASH ZT13
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/transcriptomics/data_files/liver_ZT13_gene_count_matrix.tsv
#md_file=$path/transcriptomics/liver_metadata_wnashscore_fibste.tsv
#output_path=$path/transcriptomics/rpca_results_liver/rpca_results_NASH_ZT13
#comparison=condition
#comparison2=NASH_category
#comparison3=fibrosis_stage_new
#comparison4=steatosis_grade_new

#NASH ZT1
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/transcriptomics/data_files/liver_ZT1_gene_count_matrix.tsv
#md_file=$path/transcriptomics/liver_metadata_wnashscore_fibste.tsv
#output_path=$path/transcriptomics/rpca_results_liver/rpca_results_NASH_ZT1
#comparison=condition
#comparison2=NASH_category
#comparison3=fibrosis_stage_new
#comparison4=steatosis_grade_new

#NASH FAFT (protein coding)
path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
mtb_file=$path/transcriptomics/data_files/liver_gene_count_FAFTprot_matrix.tsv
md_file=$path/transcriptomics/liver_metadata_wnashscore_fibste.tsv
output_path=$path/transcriptomics/rpca_results_liver/rpca_results_NASH_FAFTprot
comparison=condition
comparison2=timepoint
comparison3=NASH_category
comparison4=fibrosis_stage_new
comparison5=steatosis_grade_new

#all NASH (protein coding)
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/transcriptomics/data_files/liver_gene_count_prot_matrix.tsv
#md_file=$path/transcriptomics/liver_metadata_wnashscore_fibste.tsv
#output_path=$path/transcriptomics/rpca_results_liver/rpca_results_NASH_prot
#comparison=condition
#comparison2=timepoint
#comparison3=NASH_category
#comparison4=fibrosis_stage_new
#comparison5=steatosis_grade_new

#all NASH
#path=/mnt/zarrinpar/scratch/sfloresr/NASH_KF
#mtb_file=$path/transcriptomics/data_files/liver_gene_count_matrix.tsv
#md_file=$path/transcriptomics/liver_metadata_wnashscore_fibste.tsv
#output_path=$path/transcriptomics/rpca_results_liver/rpca_results_NASH_all
#comparison=condition
#comparison2=timepoint
#comparison3=NASH_category
#comparison4=fibrosis_stage_new
#comparison5=steatosis_grade_new

biom convert \
 -i $mtb_file \
 -o ${mtb_file/%.tsv/.biom} \
 --table-type="OTU table" \
 --to-hdf5

qiime tools import \
 --input-path ${mtb_file/%.tsv/.biom} \
 --output-path ${mtb_file/%.tsv/.qza} \
 --type FeatureTable[Frequency]

qiime deicode rpca \
  --i-table ${mtb_file/%.tsv/.qza}\
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
