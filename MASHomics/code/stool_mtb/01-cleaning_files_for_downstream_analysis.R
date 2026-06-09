setwd("~/Notebooks/sfloresr/MASH-TRF/MASHomics/")

library(tidyverse)
library(data.table)
library("qiime2R")

##################################################################
#inputs
mash_metadata<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/20231108_gnps_run/metadata_table/metadata_table-00000.txt"
mash_quanttab<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/20231108_gnps_run/quantification_table/quantification_table-00000.csv"
NASH_histology<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/nas_scoring_analysis_nash.csv"
mtb_annotations<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/20231108_gnps_run/DB_result/8dfc7fa4496841558d688c25f58c63c0.tsv"
mz_values<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/20231108_gnps_run/DB_result/mz_annot_v2.txt"
data_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/"

#########################################################################
#cleaning GNPS metadata files

md<-fread(mash_metadata)%>%
  filter(grepl("mzML",filename))%>%
  filter(!grepl("blank",filename))%>%
  mutate(sample_name=paste(filename,"Peak area",sep=" "),
         ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition))%>%
  dplyr::select(-c(22:29))%>%
  dplyr::select(sample_name,everything())

write.table(md,paste0(data_dir,"20231108_gnps_run/metadata_table/metadata_table-00000-clean.txt"),sep = "\t",row.names = FALSE, quote=FALSE)  

#just MASH
md_mash<-md%>%filter(ATTRIBUTE_disease_stage!="HCC")
write.table(md_mash,paste0(data_dir,"20231108_gnps_run/metadata_table/metadata_table-00000-clean-mash.txt"),sep = "\t",row.names = FALSE, quote=FALSE)  

#just HCC
md_hcc<-md%>%filter(ATTRIBUTE_disease_stage=="HCC")
write.table(md_hcc,paste0(data_dir,"20231108_gnps_run/metadata_table/metadata_table-00000-clean-hcc.txt"),sep = "\t",row.names = FALSE, quote=FALSE)  

histo <- fread(NASH_histology, header=TRUE)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         NASH_category=ifelse(NASH_category=="N/A","not_applicable",NASH_category),
         steatosis_grade=as.character(steatosis_grade),
         lobular_inflammation=as.character(lobular_inflammation),
         hepatocyte_ballooning=as.character(hepatocyte_ballooning),
         portal_inflammation=as.character(portal_inflammation),
         nucleomegallyanisonucleosis=as.character(nucleomegallyanisonucleosis),
         fibrosis_stage=as.character(fibrosis_stage),
         NASH_score=as.character(NASH_score),
         glycogenated_nuclei=factor(glycogenated_nuclei, level=c('absent', 'rare', 'few', 'many')))%>%
  mutate(NASH_category=factor(NASH_category,level=c("Non_NASH","NASH","not_applicable")),
         condition=factor(condition,level=c("NA","FA","FT")))%>%
  dplyr::select(mouse_id,NASH_category,fibrosis_stage,steatosis_grade)%>%
  dplyr::rename(ATTRIBUTE_host_subject_id=mouse_id)

md<-fread(md_mash)%>%
  mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition))%>%
  left_join(.,histo,by="ATTRIBUTE_host_subject_id")%>%
  mutate(fibrosis_stage_new=ifelse(fibrosis_stage==0,"0","1+"),
         steatosis_grade_new=ifelse(steatosis_grade==0,"0","1+"))%>%
  arrange(ATTRIBUTE_host_subject_id)

md<-md%>%
  mutate(steatosis_grade_new=ifelse((ATTRIBUTE_condition=="NA"& !is.na(NASH_category)),"none",steatosis_grade_new),
         fibrosis_stage_new=ifelse((ATTRIBUTE_condition=="NA"&!is.na(NASH_category)),"none",fibrosis_stage_new))

write.table(md,paste0(data_dir,"NASH_stool_metabolomics_metadata_new.txt"),sep = "\t",row.names = FALSE, quote=FALSE)  

#########################################################################
#cleaning GNPS quantification files

#clean & just keep mash samples
mtb<-fread(mash_quanttab)%>%
  dplyr::select(`row ID`, any_of(md_mash$sample_name))
write.table(mtb,paste0(data_dir,"20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-rmblnk.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#subsets for analysis

#FAFT wk8
mdwk8FAFT<-md%>%filter(ATTRIBUTE_host_age_at_sample_collection==7 & ATTRIBUTE_condition!="NA")
mtbw8FAFT<-mtb%>% dplyr::select(`row ID`, any_of(mdwk8FAFT$sample_name))
write.table(mtbw8FAFT,paste0("20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk8-FAFT.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FAFT wk8 ZT1
mdwk8ZT1FAFT<-md%>%filter(ATTRIBUTE_host_age_at_sample_collection==7 & ATTRIBUTE_condition!="NA" & ATTRIBUTE_timepoint=="ZT1")
mtbw8ZT1FAFT<-mtb%>% dplyr::select(`row ID`, any_of(mdwk8ZT1FAFT$sample_name))
write.table(mtbw8ZT1FAFT,paste0("20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk8-ZT1FAFT.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
mtbw8ZT1FAFT_log <- mtbw8ZT1FAFT %>% column_to_rownames("row ID")%>%decostand(method = "log")%>%
  as.data.frame()%>%rownames_to_column("row ID")
write.table(mtbw8ZT1FAFT_log, paste0("20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk8-ZT1FAFT_log.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FAFT wk8 ZT13
mdwk8ZT13FAFT<-md%>%filter(ATTRIBUTE_host_age_at_sample_collection==7 & ATTRIBUTE_condition!="NA" & ATTRIBUTE_timepoint=="ZT13")
mtbw8ZT13FAFT<-mtb%>% dplyr::select(`row ID`, any_of(mdwk8ZT13FAFT$sample_name))
write.table(mtbw8ZT13FAFT,paste0("20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk8-ZT13FAFT.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
mtbw8ZT13FAFT_log <- mtbw8ZT13FAFT %>% column_to_rownames("row ID")%>%decostand(method = "log")%>%
  as.data.frame()%>%rownames_to_column("row ID")
write.table(mtbw8ZT13FAFT_log, paste0("20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk8-ZT13FAFT_log.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FAFT wk12
mdwk12FAFT<-md%>%filter(ATTRIBUTE_host_age_at_sample_collection==12 & ATTRIBUTE_condition!="NA")
mtbw12FAFT<-mtb%>% dplyr::select(`row ID`, any_of(mdwk12FAFT$sample_name))
write.table(mtbw12FAFT,"mtb/stool/data_files/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk12-FAFT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#FAFT wk12 ZT1
mdwk12ZT1FAFT<-md%>%filter(ATTRIBUTE_host_age_at_sample_collection==12 & ATTRIBUTE_condition!="NA" & ATTRIBUTE_timepoint=="ZT1")
mtbw12ZT1FAFT<-mtb%>% dplyr::select(`row ID`, any_of(mdwk12ZT1FAFT$sample_name))
write.table(mtbw12ZT1FAFT,"mtb/stool/data_files/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk12-ZT1FAFT.txt",sep = "\t",row.names = FALSE, quote=FALSE)
mtbw12ZT1FAFT_log <- mtbw12ZT1FAFT %>% column_to_rownames("row ID")%>%decostand(method = "log")%>%
  as.data.frame()%>%rownames_to_column("row ID")
write.table(mtbw12ZT1FAFT_log, "mtb/stool/data_files/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk12-ZT1FAFT_log.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#FAFT wk12 ZT13
mdwk12ZT13FAFT<-md%>%filter(ATTRIBUTE_host_age_at_sample_collection==12 & ATTRIBUTE_condition!="NA" & ATTRIBUTE_timepoint=="ZT13")
mtbw12ZT13FAFT<-mtb%>% dplyr::select(`row ID`, any_of(mdwk12ZT13FAFT$sample_name))
write.table(mtbw12ZT13FAFT,"mtb/stool/data_files/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk12-ZT13FAFT.txt",sep = "\t",row.names = FALSE, quote=FALSE)
mtbw12ZT13FAFT_log <- mtbw12ZT13FAFT %>% column_to_rownames("row ID")%>%decostand(method = "log")%>%
  as.data.frame()%>%rownames_to_column("row ID")
write.table(mtbw12ZT13FAFT_log, "mtb/stool/data_files/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-Wk12-ZT13FAFT_log.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#########################################################################
#clean annotation table

annot<-fread(mtb_annotations)%>%
  dplyr::rename(FeatureID=`#Scan#`)
write.table(annot,paste0(data_dir,"20231108_gnps_run/DB_result/DB_results_clean.txt"),sep = "\t",row.names = FALSE, quote=FALSE)  

mz_info<-fread(mz_values)

annotations<-annot%>%
  select(FeatureID,Compound_Name)%>%
  mutate(FeatureID=as.integer(FeatureID),
         Compound_Name=gsub("\\(predic.*","",Compound_Name),
         Compound_Name=gsub('"',"",Compound_Name))%>%
  right_join(.,mz_info,by="FeatureID")
write.table(annotations,paste0(data_dir,"20231108_gnps_run/DB_result/feature_annotations_clean.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

#########################################################################
# create quantification table with extended cytoscape annotations

BA_annot<-fread(paste0(data_dir,"20231108_gnps_run/gnps_molecular_network_iin_collapse_graphml/cytoscape_BA_nodetable.csv"))%>%
  dplyr::rename(FeatureID=name)

FA_annot<-fread(paste0(data_dir,"mtb/stool/data_files/20231108_gnps_run/gnps_molecular_network_iin_collapse_graphml/cytoscape_FA_nodetable.csv"))%>%
  dplyr::rename(FeatureID=name)

AA_annot<-fread(paste0(data_dir,"20231108_gnps_run/gnps_molecular_network_iin_collapse_graphml/cytoscape_amine_nodetable.csv"))%>%
  dplyr::rename(FeatureID=name)

gly_annot<-fread(paste0(data_dir,"20231108_gnps_run/gnps_molecular_network_iin_collapse_graphml/cytoscape_glycol_nodetable.csv"))%>%
  dplyr::rename(FeatureID=name)

car_annot<-fread(paste0(data_dir,"20231108_gnps_run/gnps_molecular_network_iin_collapse_graphml/cytoscape_carnitine_nodetable.csv"))%>%
  dplyr::rename(FeatureID=name)

phos_annot<-fread(paste0(data_dir,"20231108_gnps_run/gnps_molecular_network_iin_collapse_graphml/cytoscape_phosphocholine_nodetable.csv"))%>%
  dplyr::rename(FeatureID=name)

bil_annot<-fread(paste0(data_dir,"20231108_gnps_run/gnps_molecular_network_iin_collapse_graphml/cytoscape_bilirubin_nodetable.csv"))%>%
  dplyr::rename(FeatureID=name)

mz<-fread(paste0(data_dir,"20231108_gnps_run/quantification_table/quantification_table-00000.csv"))%>%
  dplyr::select(c(1:3))%>%
  dplyr::rename(FeatureID=`row ID`,mz= `row m/z`,RT=`row retention time`)%>%
  mutate(mtb_group = dplyr::case_when(FeatureID %in% BA_annot$FeatureID ~ "bile acids", 
                                      FeatureID %in% FA_annot$FeatureID ~ "fatty acids & other lipids",
                                      FeatureID %in% AA_annot$FeatureID ~ "amines",
                                      FeatureID %in% gly_annot$FeatureID ~ "glycols",
                                      FeatureID %in% car_annot$FeatureID ~ "carnitines",
                                      FeatureID %in% phos_annot$FeatureID ~ "phosphocholines",
                                      FeatureID %in% bil_annot$FeatureID ~ "bilirubins",
                                      TRUE ~ "unknown"))
write.table(mz,paste0(data_dir,"20231108_gnps_run/DB_result/mz_annot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)  
write.table(mz,paste0(data_dir,"20231108_gnps_run/DB_result/mz_annot_v2.txt"),sep = "\t",row.names = FALSE, quote=FALSE)  

mtb<-fread(paste0(data_dir,"20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-rmblnk.txt"))

#make long tabl
mtb_wmd<-mtb%>%dplyr::rename(FeatureID=`row ID`)%>%
  gather(sample_name,counts,-FeatureID)%>%
  left_join(.,mz,by="FeatureID")%>%
  left_join(.,md,by="sample_name")%>%
  filter(!is.na(NASH_category))%>%
  left_join(.,annot,by="FeatureID")

#write.table(mtb_wmd,paste0(data_dir,"NASH_stool_metabolomics_quanttab_long_wmd.txt"),sep = "\t",row.names = FALSE, quote=FALSE)  
write.table(mtb_wmd,paste0(data_dir,"NASH_stool_metabolomics_quanttab_long_wmd_v2.txt"),sep = "\t",row.names = FALSE, quote=FALSE)  

#micromasst was run separatelly
