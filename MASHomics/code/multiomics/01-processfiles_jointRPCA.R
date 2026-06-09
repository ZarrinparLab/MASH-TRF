setwd("~/Notebooks/sfloresr/MASH-TRF/MASHomics/")

library(tidyverse)
library(data.table)
library("qiime2R")
library(VennDiagram)
library(gplots)

####################################################################
#inputs
md_16s<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_16s/s16s_metadata_cln_addmoreNASHcat.txt"
annot_16s<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_16s/mashstool16s_preprocessed_20211020_ID_13785_gg2/taxonomy.tsv"
md_mtb<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/NASH_stool_metabolomics_metadata_new.txt"
annot_mtb<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/20231108_gnps_run/DB_result/feature_annotations_clean.txt"
feattab_16s_wk12<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_16s/mashstool16s_preprocessed_20211020_ID_13785_gg2/NASH_NASH_taxonomy_filtered.gg2.asv.counts.qza"
feattab_16s_wk8<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_16s/mashstool16s_preprocessed_20211020_ID_13785_gg2/NASH_TRF_taxonomy_filtered.gg2.asv.counts.qza"
feattab_mtb<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-rmblnk.txt"
bdm_16s_res<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/stool_16s/birdman_results/"
bdm_mtb_res<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/stool_mtb/birdman_results/"
data_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/multiomics/"
####################################################################
#cleaning metadata

#clean up 16S md file (just HFD Wk12)
md_16s_wk12<-fread(md_16s)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="NA" & collection_disease_stage=="NASH")%>%
  mutate(sample_id_mouse=gsub('2$','13',sample_name))%>%
  filter(!(sample_id_mouse %in% c("18R.N1", "13L.N13","13N.N13","14R.N13","17R.N13","19R.N13")))

set.seed(1234)
data_mod <- md_16s_wk12 %>% group_by(NASH_category) %>% sample_frac(.3)

md_16s_wk12<-md_16s_wk12%>%
  mutate(traintest=ifelse(sample_id_mouse %in% data_mod$sample_id_mouse, "test","train"))%>%
  arrange(sample_id_mouse)

write.table(md_16s_wk12, paste0(data_dir,"jrpca_16s_key.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#clean up mtb md file (just HFD Wk12)

md_mtb_wk12<-fread(md_mtb)%>%
  mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition))%>%
  filter(ATTRIBUTE_condition!="NA" & ATTRIBUTE_disease_stage=="NASH")%>%
  filter(!is.na(NASH_category))%>%
  mutate(traintest=ifelse(sample_id %in% data_mod$sample_id_mouse, "test","train"))%>%
  arrange(sample_id)

# list_venn <- list(md_16s$sample_id_mouse,md_mtb$sample_id_mouse)
# ItemsList <- venn(list_venn, show.plot = FALSE)
# all<-attributes(ItemsList)$intersections

write.table(md_mtb_wk12, paste0(data_dir,"jrpca_mtb_key.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#metadata for both 16s and mtb (wk12)
jrpca_md<-md_16s_wk12%>%
  dplyr::select(sample_id_mouse,host_subject_id,condition,collection_timepoint,collection_disease_stage,
                NASH_category,steatosis_grade_new,fibrosis_stage_new,traintest)%>%
  dplyr::rename(sample_name=sample_id_mouse)%>%
  mutate(cond_NASH=paste(condition,NASH_category, sep="_"))

write.table(jrpca_md, paste0(data_dir,"jrpca_metadata.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#clean up 16S md file (just HFD Wk8)
md_16s_wk8<-fread(md_16s)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="NA" & collection_disease_stage=="TRF")%>%
  mutate(sample_id_mouse=gsub('2$','13',sample_name))%>%
  filter(!(sample_id_mouse %in% c("18R.T1","19R.T1","19R.T13")))

set.seed(1234)
data_mod <- md_16s_wk8 %>% group_by(NASH_category) %>% sample_frac(.3)

md_16s_wk8<-md_16s_wk8%>%
  mutate(traintest=ifelse(sample_id_mouse %in% data_mod$sample_id_mouse, "test","train"))%>%
  arrange(sample_id_mouse)

write.table(md_16s_wk8, paste0(data_dir,"jrpca_16s_key_wk8.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#clean up mtb md file (just HFD Wk8)
md_mtb_wk8<-fread(md_mtb)%>%
  mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition))%>%
  filter(ATTRIBUTE_condition!="NA" & ATTRIBUTE_disease_stage=="TRF")%>%
  filter(!is.na(NASH_category))%>%
  dplyr::rename(sample_id_mouse=sample_id)%>%
  mutate(traintest=ifelse(sample_id_mouse %in% data_mod$sample_id_mouse, "test","train"))%>%
  arrange(sample_id_mouse)
write.table(md_mtb_wk8, paste0(data_dir,"jrpca_mtb_key_wk8.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

# list_venn <- list(md_16s$sample_id_mouse,md_mtb$sample_id_mouse)
# ItemsList <- venn(list_venn, show.plot = FALSE)
# all<-attributes(ItemsList)$intersections

#metadata for both 16s and mtb (wk 8)
jrpca_md<-md_16s_wk8%>%
  dplyr::select(sample_id_mouse,host_subject_id,condition,collection_timepoint,collection_disease_stage,
                NASH_category,steatosis_grade_new,fibrosis_stage_new,traintest)%>%
  dplyr::rename(sample_name=sample_id_mouse)%>%
  mutate(cond_NASH=paste(condition,NASH_category, sep="_"))

write.table(jrpca_md, paste(data_dir,"jrpca_metadata_wk8.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

####################################################################
#clean annotations

#make 16s feature/annotation key
m16s_annot<-fread(annot_16s)%>%
  separate(Taxon, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep=";", remove=FALSE)

write.table(m16s_annot, paste0(data_dir,"jrpca_16s_annotation_key_gg2.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#make mtb feature/annotation key
class<-fread(paste0(data_dir,"canopus_summary.tsv"))%>%
  separate(name, c(NA,"FeatureID"),sep="sirius_")%>%
  mutate(FeatureID=as.integer(FeatureID))
mtb_annot<-fread(annot_mtb)%>%
  left_join(.,class,by="FeatureID")

write.table(mtb_annot, paste0(data_dir,"jrpca_mtb_annotation_key_wcanopus.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

####################################################################
#clean feature tables

#clean up 16S feattab file (just HFD Wk12)
m16s<-read_qza(feattab_16s_wk12)$data%>%
  as.data.frame()%>%rownames_to_column("FeatureID")%>%
  dplyr::select(FeatureID,all_of(md_16s_wk12$sample_name))

names(m16s) <- c("FeatureID","13L.N1","13N.N1","14L.N1","14L.N13","14R.N1","15L.N1","15L.N13","16R.N1",   
                 "16R.N13","17N.N1","17N.N13","17R.N1","18L.N1","18L.N13","18R.N13","19R.N1","20N.N1","20N.N13")
write.table(m16s, paste0(data_dir,"jrpca_16s_feattab_gg2.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#clean up mtb feattab file (just HFD Wk12)
mtb<-fread(feattab_mtb)%>%
  as.data.frame()%>%dplyr::rename(FeatureID=`row ID`)%>%
  dplyr::select(FeatureID,all_of(md_mtb_wk12$sample_name))

names(mtb) <- c("FeatureID","13L.N1","13N.N1","14L.N1","14L.N13","14R.N1","15L.N1","15L.N13","16R.N1",   
                "16R.N13","17N.N1","17N.N13","17R.N1","18L.N1","18L.N13","18R.N13","19R.N1","20N.N1","20N.N13")
write.table(mtb, paste0(data_dir,"jrpca_mtb_feattab.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#clean up 16S feattab file (just HFD Wk8)
m16s<-read_qza(feattab_16s_wk8)$data%>%
  as.data.frame()%>%rownames_to_column("FeatureID")%>%
  dplyr::select(FeatureID,all_of(md_16s_wk8$sample_name))

names(m16s) <- c("FeatureID","13L.T1","13L.T13","13N.T1","13N.T13","14L.T1","14L.T13","14R.T1","14R.T13","15L.T1",  
                 "15L.T13","16R.T1","16R.T13","17N.T1","17N.T13","17R.T1","17R.T13","18L.T1","18L.T13","18R.T13",
                 "20N.T1","20N.T13")
write.table(m16s, paste0(data_dir,"jrpca_16s_feattab_wk8_gg2.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#clean up mtb feattab file (just HFD Wk12)
mtb<-fread(feattab_mtb)%>%
  as.data.frame()%>%dplyr::rename(FeatureID=`row ID`)%>%
  dplyr::select(FeatureID,all_of(md_mtb_wk8$sample_name))

names(mtb) <- c("FeatureID","13L.T1","13L.T13","13N.T1","13N.T13","14L.T1","14L.T13","14R.T1","14R.T13","15L.T1","15L.T13",
                "16R.T1","16R.T13","17N.T1","17N.T13","17R.T1","17R.T13","18L.T1","18L.T13","18R.T13","20N.T1","20N.T13")
write.table(mtb, paste0(data_dir,"jrpca_mtb_feattab_wk8.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

####################################################################

#create a table with just the microbes (16s) that were different using birdman

FA12<-fread(paste0(bdm_16s_res,"birdman_justcred_wk12FA_nonNASHNASH_results_gg2.tsv"))
FT12<-fread(paste0(bdm_16s_res,"birdman_justcred_wk12FT_nonNASHNASH_results_gg2.tsv"))

list_venn <- list(FA12$FeatureID,FT12$FeatureID)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-unlist(attributes(ItemsList)$intersections)

m16s_bdm<-fread(paste0(data_dir,"jrpca_16s_feattab.txt"))%>%
  filter(FeatureID %in% all) #47 microbes

write.table(m16s_bdm, paste0(data_dir,"jrpca_16sbdm_feattab_gg2.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#create a table of microbes that are present in week 8 and 12 

FA8<-fread(paste0(bdm_16s_res,"birdman_justcred_wk8FA_nonNASHNASH_results_gg2.tsv"))
FT8<-fread(paste0(bdm_16s_res,"birdman_justcred_wk8FT_nonNASHNASH_results_gg2.tsv"))

list_venn <- list(FA8$FeatureID,FA12$FeatureID)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

m6FA<-all$`A:B`

list_venn2 <- list(FT8$FeatureID,FT12$FeatureID)
ItemsList2 <- venn(list_venn2, show.plot = FALSE)
all2<-attributes(ItemsList2)$intersections

m3FT<-all2$`A:B`

list_venn3 <- list(m6FA,m3FT)
ItemsList3 <- venn(list_venn3, show.plot = FALSE)
all3<-unlist(attributes(ItemsList3)$intersections)

#W12
m16s_shared_bdm<-fread(paste0(data_dir,"jrpca_16s_feattab_gg2.txt"))%>%
  filter(FeatureID %in% all3) #10 microbes

write.table(m16s_shared_bdm, paste0(data_dir,"jrpca_16sbdm_shared_feattab_gg2.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#W8
m16s_shared_bdm_w8<-fread(paste0(data_dir,"jrpca_16s_feattab_wk8_gg2.txt"))%>%
  filter(FeatureID %in% all3) #10 microbes

write.table(m16s_shared_bdm_w8, "multiomics_16smtb/jrpca_16sbdm_shared_feattab_wk8_gg2.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#create a table with just the mtbs that where diff expr from bdm (in both wk8 and wk12)
wk8ZT1<-fread(paste0(bdm_mtb_res,"birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv"))
wk12ZT1<-fread(paste0(bdm_mtb_res,"birdman_justcred_wk12ZT1_nonNASHNASH_results.tsv"))
list_venn <- list(Wk8 = wk8ZT1$Name,
                  Wk12 = wk12ZT1$Name)
ItemsList <- venn(list_venn, show.plot = FALSE)
all1<-attributes(ItemsList)$intersections

wk8ZT13<-fread(paste0(bdm_mtb_res,"birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv"))
wk12ZT13<-fread(paste0(bdm_mtb_res,"birdman_justcred_wk12ZT13_nonNASHNASH_results.tsv"))

list_venn <- list(Wk8 = wk8ZT13$Name,
                  Wk12 = wk12ZT13$Name)
ItemsList <- venn(list_venn, show.plot = FALSE)
all2<-attributes(ItemsList)$intersections

list_venn <- list(ZT1 = all1$`Wk8:Wk12`,
                  ZT13 = all2$`Wk8:Wk12`)
ItemsList4 <- venn(list_venn, show.plot = FALSE)
all4<-attributes(ItemsList4)$intersections

mtb_list<-c(all1$`Wk8:Wk12`,all2$`Wk8:Wk12`)%>%unique()  #41 unique 
mtb_list <- sapply(strsplit(mtb_list, " "), function(x) x[1])

mtb_sel41<-fread(paste0(data_dir,"jrpca_mtb_feattab.txt"))%>% 
  filter(FeatureID %in% mtb_list) 

write.table(mtb_sel41, paste0(data_dir,"jrpca_mtbbdm_sel41_feattab.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

mtb_sel41_w8<-fread(paste0(data_dir,"jrpca_mtb_feattab_wk8.txt"))%>%
  filter(FeatureID %in% mtb_list) 

write.table(mtb_sel41_w8, paste0(data_dir,"jrpca_mtbbdm_sel41_feattab_wk8.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

####################################################################
#create 16s and mtb keys for jrpca plotting

sel16s_NASHcat<-data.frame(FeatureID=all3,NASH_bdmcat=c("Non-MASH","Non-MASH","Non-MASH","Non-MASH",
                                                        "Non-MASH:MASH","MASH","Non-MASH","Non-MASH","Non-MASH:MASH","Non-MASH:MASH"),
                           condition=c("FA","FA","FA","FA","FA","FA","FT","FT","FT","FA:FT"))

write.table(sel16s_NASHcat, paste0(data_dir,"jrpca_sel16s_NASHcat.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

sel41_ZTcat<-data.frame(name=c(all4$ZT1,all4$ZT13,all4$`ZT1:ZT13`),
                        ZT_cat=c(rep("ZT1",26),rep("ZT13",12),rep("ZT1:ZT13",3)))%>%
  mutate(FeatureID = gsub(" .*", "", name))%>%
  dplyr::select(FeatureID,ZT_cat)

write.table(sel41_ZTcat, paste0(data_dir,"jrpca_sel41_ZTcat.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
