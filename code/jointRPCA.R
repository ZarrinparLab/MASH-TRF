setwd("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/")

library(tidyverse)
library(data.table)
library("qiime2R")
library("Biostrings")
library(ggcorrplot)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(scatterplot3d)
library(VennDiagram)
library(gplots)
####################################################################

#clean up 16S md file (just HFD Wk12)

md_16s<-fread("16s/stool/metadata_cln_addmoreNASHcat.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="NA" & collection_disease_stage=="NASH")%>%
  mutate(sample_id_mouse=gsub('2$','13',sample_name))%>%
  filter(!(sample_id_mouse %in% c("18R.N1", "13L.N13","13N.N13","14R.N13","17R.N13","19R.N13")))%>%
  mutate(traintest=ifelse(sample_id_mouse %in% data_mod$sample_id_mouse, "test","train"))%>%
  arrange(sample_id_mouse)

write.table(md_16s, "multiomics_16smtb/jrpca_16s_key.txt",sep = "\t",row.names = FALSE, quote=FALSE)

set.seed(1234)
data_mod <- md_16s %>% group_by(NASH_category) %>% sample_frac(.3)

#clean up mtb md file (just HFD Wk12)
md_mtb<-fread("mtb/stool/NASH_stool_metabolomics_metadata_new.txt")%>%
  mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition))%>%
  filter(ATTRIBUTE_condition!="NA" & ATTRIBUTE_disease_stage=="NASH")%>%
  filter(!is.na(NASH_category))%>%
  mutate(traintest=ifelse(sample_id %in% data_mod$sample_id_mouse, "test","train"))%>%
  arrange(sample_id)

# list_venn <- list(md_16s$sample_id_mouse,md_mtb$sample_id_mouse)
# ItemsList <- venn(list_venn, show.plot = FALSE)
# all<-attributes(ItemsList)$intersections

write.table(md_mtb, "multiomics_16smtb/jrpca_mtb_key.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#metadata for both 16s and mtb (wk12)
jrpca_md<-md_16s%>%
  dplyr::select(sample_id_mouse,host_subject_id,condition,collection_timepoint,collection_disease_stage,
                NASH_category,steatosis_grade_new,fibrosis_stage_new,traintest)%>%
  dplyr::rename(sample_name=sample_id_mouse)%>%
  mutate(cond_NASH=paste(condition,NASH_category, sep="_"))

write.table(jrpca_md, "multiomics_16smtb/jrpca_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#clean up 16S feattab file (just HFD Wk12)
m16s<-read_qza("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/NASH_NASH_taxonomy_filtered.gg2.asv.counts.qza")$data%>%
  as.data.frame()%>%rownames_to_column("FeatureID")%>%
  dplyr::select(FeatureID,all_of(md_16s$sample_name))

names(m16s) <- c("FeatureID","13L.N1","13N.N1","14L.N1","14L.N13","14R.N1","15L.N1","15L.N13","16R.N1",   
                 "16R.N13","17N.N1","17N.N13","17R.N1","18L.N1","18L.N13","18R.N13","19R.N1","20N.N1","20N.N13")
write.table(m16s, "multiomics_16smtb/jrpca_16s_feattab_gg2.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#make 16s feature/annotation key
m16s_annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/taxonomy.tsv")%>%
  separate(Taxon, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep=";", remove=FALSE)

write.table(m16s_annot, "multiomics_16smtb/jrpca_16s_annotation_key_gg2.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#clean up mtb feattab file (just HFD Wk12)
#mtb<-read_qza("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/metabolomics/stool/data/filtered_featuretable_allzt_NASH_NASH.qza")$data%>%
mtb<-fread("mtb/stool/data_files/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-rmblnk.txt")%>%
  as.data.frame()%>%dplyr::rename(FeatureID=`row ID`)%>%
  dplyr::select(FeatureID,all_of(md_mtb$sample_name))

names(mtb) <- c("FeatureID","13L.N1","13N.N1","14L.N1","14L.N13","14R.N1","15L.N1","15L.N13","16R.N1",   
                "16R.N13","17N.N1","17N.N13","17R.N1","18L.N1","18L.N13","18R.N13","19R.N1","20N.N1","20N.N13")
write.table(mtb, "multiomics_16smtb/jrpca_mtb_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#make mtb feature/annotation key
class<-fread("mtb/stool/old_files/canopus_summary.tsv")%>%
  separate(name, c(NA,"FeatureID"),sep="sirius_")%>%
  mutate(FeatureID=as.integer(FeatureID))
mtb_annot<-fread("mtb/stool/data_files/20231108_gnps_run/DB_result/feature_annotations_clean.txt")%>%
  left_join(.,class,by="FeatureID")

write.table(mtb_annot, "multiomics_16smtb/jrpca_mtb_annotation_key_wcanopus.txt",sep = "\t",row.names = FALSE, quote=FALSE)

####################################################################

#clean up 16S md file (just HFD Wk8)
md_16s<-fread("16s/stool/metadata_cln_addmoreNASHcat.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="NA" & collection_disease_stage=="TRF")%>%
  mutate(sample_id_mouse=gsub('2$','13',sample_name))%>%
  filter(!(sample_id_mouse %in% c("18R.T1","19R.T1","19R.T13")))%>%
  mutate(traintest=ifelse(sample_id_mouse %in% data_mod$sample_id_mouse, "test","train"))%>%
  arrange(sample_id_mouse)

write.table(md_16s, "multiomics_16smtb/jrpca_16s_key_wk8.txt",sep = "\t",row.names = FALSE, quote=FALSE)

set.seed(1234)
data_mod <- md_16s %>% group_by(NASH_category) %>% sample_frac(.3)

#clean up mtb md file (just HFD Wk8)
md_mtb<-fread("mtb/stool/NASH_stool_metabolomics_metadata_new.txt")%>%
  mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition))%>%
  filter(ATTRIBUTE_condition!="NA" & ATTRIBUTE_disease_stage=="TRF")%>%
  filter(!is.na(NASH_category))%>%
  dplyr::rename(sample_id_mouse=sample_id)%>%
  mutate(traintest=ifelse(sample_id_mouse %in% data_mod$sample_id_mouse, "test","train"))%>%
  arrange(sample_id_mouse)
write.table(md_mtb, "multiomics_16smtb/jrpca_mtb_key_wk8.txt",sep = "\t",row.names = FALSE, quote=FALSE)

# list_venn <- list(md_16s$sample_id_mouse,md_mtb$sample_id_mouse)
# ItemsList <- venn(list_venn, show.plot = FALSE)
# all<-attributes(ItemsList)$intersections


#metadata for both 16s and mtb
jrpca_md<-md_16s%>%
  dplyr::select(sample_id_mouse,host_subject_id,condition,collection_timepoint,collection_disease_stage,
                NASH_category,steatosis_grade_new,fibrosis_stage_new,traintest)%>%
  dplyr::rename(sample_name=sample_id_mouse)%>%
  mutate(cond_NASH=paste(condition,NASH_category, sep="_"))

write.table(jrpca_md, "multiomics_16smtb/jrpca_metadata_wk8.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#clean up 16S feattab file (just HFD Wk8)
m16s<-read_qza("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/NASH_TRF_taxonomy_filtered.gg2.asv.counts.qza")$data%>%
  as.data.frame()%>%rownames_to_column("FeatureID")%>%
  dplyr::select(FeatureID,all_of(md_16s$sample_name))

names(m16s) <- c("FeatureID","13L.T1","13L.T13","13N.T1","13N.T13","14L.T1","14L.T13","14R.T1","14R.T13","15L.T1",  
                 "15L.T13","16R.T1","16R.T13","17N.T1","17N.T13","17R.T1","17R.T13","18L.T1","18L.T13","18R.T13",
                 "20N.T1","20N.T13")
write.table(m16s, "multiomics_16smtb/jrpca_16s_feattab_wk8_gg2.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#clean up mtb feattab file (just HFD Wk12)
mtb<-fread("mtb/stool/data_files/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-rmblnk.txt")%>%
  as.data.frame()%>%dplyr::rename(FeatureID=`row ID`)%>%
  dplyr::select(FeatureID,all_of(md_mtb$sample_name))

names(mtb) <- c("FeatureID","13L.T1","13L.T13","13N.T1","13N.T13","14L.T1","14L.T13","14R.T1","14R.T13","15L.T1","15L.T13",
                "16R.T1","16R.T13","17N.T1","17N.T13","17R.T1","17R.T13","18L.T1","18L.T13","18R.T13","20N.T1","20N.T13")
write.table(mtb, "multiomics_16smtb/jrpca_mtb_feattab_wk8.txt",sep = "\t",row.names = FALSE, quote=FALSE)


####################################################################

#create a table with just the microbes (16s) that were different using birdman

FA12<-fread("16s/stool/birdman_results/birdman_justcred_wk12FA_nonNASHNASH_results_gg2.tsv")
FT12<-fread("16s/stool/birdman_results/birdman_justcred_wk12FT_nonNASHNASH_results_gg2.tsv")

list_venn <- list(FA12$FeatureID,FT12$FeatureID)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-unlist(attributes(ItemsList)$intersections)

m16s_bdm<-fread("multiomics_16smtb/jrpca_16s_feattab.txt")%>%
  filter(FeatureID %in% all) #47 microbes

write.table(m16s_bdm, "multiomics_16smtb/jrpca_16sbdm_feattab_gg2.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#create a table of microbes that are present in week 8 and 12 

FA8<-fread("16s/stool/birdman_results/birdman_justcred_wk8FA_nonNASHNASH_results_gg2.tsv")
FT8<-fread("16s/stool/birdman_results/birdman_justcred_wk8FT_nonNASHNASH_results_gg2.tsv")

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
m16s_shared_bdm<-fread("multiomics_16smtb/jrpca_16s_feattab_gg2.txt")%>%
  filter(FeatureID %in% all3) #10 microbes

write.table(m16s_shared_bdm, "multiomics_16smtb/jrpca_16sbdm_shared_feattab_gg2.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#W8
m16s_shared_bdm_w8<-fread("multiomics_16smtb/jrpca_16s_feattab_wk8_gg2.txt")%>%
  filter(FeatureID %in% all3) #10 microbes

write.table(m16s_shared_bdm_w8, "multiomics_16smtb/jrpca_16sbdm_shared_feattab_wk8_gg2.txt",sep = "\t",row.names = FALSE, quote=FALSE)


#create a table with just the mtbs that mapped to bile acids 
class<-fread("mtb/stool/data_files/20231108_gnps_run/DB_result/feature_annotations_clean.txt")
BA<-class%>%filter(mtb_group=="bile acids")
mtb_BA<-fread("multiomics_16smtb/jrpca_mtb_feattab.txt")%>%
  filter(FeatureID %in% BA$FeatureID) #290 mtbs

write.table(mtb_BA, "multiomics_16smtb/jrpca_mtbBA_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)


#create a table with just the mtbs that where diff expr from bdm (in both wk8 and wk12)
wk8ZT1<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv")
wk12ZT1<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk12ZT1_nonNASHNASH_results.tsv")
list_venn <- list(Wk8 = wk8ZT1$Name,
                  Wk12 = wk12ZT1$Name)
ItemsList <- venn(list_venn, show.plot = FALSE)
all1<-attributes(ItemsList)$intersections

wk8ZT13<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv")
wk12ZT13<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk12ZT13_nonNASHNASH_results.tsv")

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

mtb_sel41<-fread("multiomics_16smtb/jrpca_mtb_feattab.txt")%>%
  filter(FeatureID %in% mtb_list) 

write.table(mtb_sel41, "multiomics_16smtb/jrpca_mtbbdm_sel41_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtb_sel41_w8<-fread("multiomics_16smtb/jrpca_mtb_feattab_wk8.txt")%>%
  filter(FeatureID %in% mtb_list) 

write.table(mtb_sel41_w8, "multiomics_16smtb/jrpca_mtbbdm_sel41_feattab_wk8.txt",sep = "\t",row.names = FALSE, quote=FALSE)
####################################################################
#create & plot the correlation (wk12)

sel16s_NASHcat<-data.frame(FeatureID=all3,NASH_bdmcat=c("Non-MASH","Non-MASH","Non-MASH","Non-MASH",
                                                        "Non-MASH:MASH","MASH","Non-MASH","Non-MASH","Non-MASH:MASH","Non-MASH:MASH"),
                           condition=c("FA","FA","FA","FA","FA","FA","FT","FT","FT","FA:FT"))

sel41_ZTcat<-data.frame(name=c(all4$ZT1,all4$ZT13,all4$`ZT1:ZT13`),
                        ZT_cat=c(rep("ZT1",26),rep("ZT13",12),rep("ZT1:ZT13",3)))%>%
  mutate(FeatureID = gsub(" .*", "", name))%>%
  dplyr::select(FeatureID,ZT_cat)

m16s_key<-fread("multiomics_16smtb/jrpca_16s_annotation_key_gg2.txt")%>%
  right_join(.,sel16s_NASHcat,by="FeatureID")%>%
  mutate(name=paste(FeatureID,Taxon))%>%
  dplyr::select(FeatureID,name)
mtb_key<-fread("multiomics_16smtb/jrpca_mtb_annotation_key_wcanopus.txt")%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  right_join(.,sel41_ZTcat,by="FeatureID")%>%
  mutate(name=ifelse(is.na(Compound_Name),
                     paste("mtb",FeatureID,"m/z",signif(mz,4), sep=" "),
                     paste("mtb",FeatureID,"m/z",signif(mz,4), Compound_Name,sep=" ")),
         FeatureID=as.character(FeatureID))%>%
  dplyr::select(FeatureID,name)

key<-rbind(m16s_key,mtb_key,fill=TRUE)

corr_n<-fread("multiomics_16smtb/result_bdm_sel41mtb_sel10micro_gg2/correlation_table/Correlation.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=featureid)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(53,2:11))%>%
  column_to_rownames("name")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(53,12:52))%>%
  column_to_rownames("name")

write.table(corr_n, "multiomics_16smtb/result_bdm_sel41mtb_sel10micro_gg2/correlation_table/Correlation_clean.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

p<-pheatmap(corr_n,color= inferno(10))
ggsave("multiomics_16smtb/result_bdm_sel41mtb_sel10micro_gg2/correlation_table/SFR24_0212_heatmap_annot.pdf", p,width = 20, height = 26)

m16s_key_ext<-fread("multiomics_16smtb/jrpca_16s_annotation_key_gg2.txt")%>%
  right_join(.,sel16s_NASHcat,by="FeatureID")%>%
  mutate(name=paste(FeatureID,Taxon),
         NASH_bdmcat=factor(NASH_bdmcat,levels=c("Non-MASH","MASH","Non-MASH:MASH")),
         condition=factor(condition,levels=c("FA","FT","FA:FT")))%>%
  dplyr::select(name,NASH_bdmcat,condition)%>%
  column_to_rownames("name")
mtb_key_ext<-fread("multiomics_16smtb/jrpca_mtb_annotation_key.txt")%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  right_join(.,sel41_ZTcat,by="FeatureID")%>%
  mutate(name=ifelse(is.na(Compound_Name),
                     paste("mtb",FeatureID,"m/z",signif(mz,4), sep=" "),
                     paste("mtb",FeatureID,"m/z",signif(mz,4), Compound_Name,sep=" ")),
         FeatureID=as.character(FeatureID),
         mtb_group=factor(mtb_group,levels=c("bile acids","fatty acids & other lipids","unknown")),
         ZT_cat=factor(ZT_cat,levels=c("ZT1","ZT13","ZT1:ZT13")))%>%
  dplyr::select(name,mtb_group,ZT_cat)%>%
  column_to_rownames("name")

ann_colors <- list(
  NASH_bdmcat = c("Non-MASH" = "#0000a7", "MASH" = "#c1272d", "Non-MASH:MASH" = "#92278F"),
  condition = c("FA" = "#D55E00", "FT" = "#009E73", "FA:FT" = "grey50"),
  mtb_group = c("bile acids" = "#006838", "fatty acids & other lipids" = "#F15A29", "unknown" = "#808285"),
  ZT_cat = c("ZT1" = "grey90", "ZT13" = "grey10", "ZT1:ZT13" = "grey50")
)

p<-pheatmap(corr_n, annotation_row = m16s_key_ext,annotation_col = mtb_key_ext,
            annotation_colors = ann_colors,color= inferno(10), main="Wk12")
ggsave("multiomics_16smtb/result_bdm_sel41mtb_sel10micro_gg2/correlation_table/SFR24_0212_heatmap_wbdmannot.pdf", p,width = 23, height = 26)

ordw12_micro<-rownames(corr_n[p$tree_row[["order"]],]) #microbe order
cordw12_mtb<-colnames(corr_n[,p$tree_col[["order"]]]) #mtb order
#######################################################################################

#create & plot the correlation (wk8)

corr_n<-fread("multiomics_16smtb/result_bdm_sel41mtb_sel10micro_wk8_gg2/correlation_table/Correlation.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=featureid)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(53,2:11))%>%
  column_to_rownames("name")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(53,12:52))%>%
  column_to_rownames("name")

write.table(corr_n, "multiomics_16smtb/result_bdm_sel41mtb_sel10micro_wk8_gg2/correlation_table/Correlation_clean.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

p<-pheatmap(corr_n,color= inferno(10))
ggsave("multiomics_16smtb/result_bdm_sel41mtb_sel10micro_wk8_gg2/correlation_table/SFR24_0207_heatmap_annot.pdf", p,width = 20, height = 26)

p<-pheatmap(corr_n, annotation_row = m16s_key_ext,annotation_col = mtb_key_ext,
            annotation_colors = ann_colors,color= inferno(10), main="Wk8")
ggsave("multiomics_16smtb/result_bdm_sel41mtb_sel10micro_wk8_gg2/correlation_table/SFR24_0207_heatmap_wbdmannot.pdf", p,width = 23, height = 26)


#keep the order of wk12 heatmap (turn off clustering)
corr_n<-corr_n%>%dplyr::select(all_of(cordw12_mtb))%>%
  rownames_to_column("micro")%>%
  mutate(micro =  factor(micro, levels = ordw12_micro)) %>%
  arrange(micro)%>%
  column_to_rownames("micro")
p<-pheatmap(corr_n, annotation_row = m16s_key_ext,annotation_col = mtb_key_ext, cluster_rows=FALSE, cluster_cols=FALSE,
            annotation_colors = ann_colors,color= inferno(10), main="Wk8")
ggsave("multiomics_16smtb/result_bdm_sel41mtb_sel10micro_wk8_gg2/correlation_table/SFR24_0214_heatmap_wbdmannot_noclust.pdf", p,width = 23, height = 25.5)


