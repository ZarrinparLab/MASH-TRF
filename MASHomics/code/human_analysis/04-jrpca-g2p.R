###################################################################
#create input file for joint-rpca

#16S
##metadata
m16S_md<-fread("NAFLD_human_val/caussy_nafld_16s/NAFLD_metadata.tsv")%>%filter(groups=="g1p" | groups=="g2p")%>%
  mutate(sample_name=gsub("^11635\\.", "",sample_name))%>%
  mutate(traintest=ifelse(sample_name %in% data_mod$sample_name, "test","train"))
write.table(md,"NAFLD_human_val/caussy_nafld_16s/NAFLD_metadata_clean_g2p.tsv",sep = "\t",row.names = FALSE,quote=FALSE)

set.seed(1234)
data_mod <- m16S_md %>% group_by(groups) %>% sample_frac(.3)

##quant_tab
m16S<-fread("NAFLD_human_val/caussy_nafld_16s/04.taxonomy_filtered.asv.counts.g1pg2p.txt")%>%
  select(FeatureID,all_of(m16S_md$sample_name))
write.table(m16S,"NAFLD_human_val/caussy_nafld_16s/m16S_counts_g1pg2p_jrpca.txt",sep = "\t",row.names = FALSE,quote=FALSE)

lachno_bdm<-fread("NAFLD_human_val/caussy_nafld_16s/birdman/birdman_justcred_g1pg2p_results_cred1.tsv")%>%
  filter(grepl("f__Lachnospiraceae", Name))
m16s_sub<-m16S%>%filter(FeatureID %in% lachno_bdm$FeatureID)
write.table(m16s_sub,"NAFLD_human_val/caussy_nafld_16s/m16S_counts_g1pg2p_jrpca_lachnohits.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#mtb
##metadata
mtb_md<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/metadata/clean_fecal_metadata.tsv")%>%
  mutate(sampleid=gsub("\\-","\\.",sampleid))%>%
  filter(ATTRIBUTE_groups=="G1P" | ATTRIBUTE_groups=="G2P")%>%
  filter(sampleid %in% m16S_md$sample_name)%>%
  filter(!(filename %in% c("TWDS002_RB6_01_29701.mzXML","TWCP-001_RC11_01_29719.mzXML")))

# setdiff(m16S_md$sample_name, mtb_md$sampleid)
# setdiff(mtb_md$sampleid, m16S_md$sample_name)

##quant_tab
mtb<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/quantification_table/quantification_table-G1PG2P.txt")
col_map <- setNames(paste0(mtb_md$filename, " Peak area"), mtb_md$sampleid)
mtb <- mtb %>% select(`row ID`, all_of(col_map))
write.table(mtb,"NAFLD_human_val/caussy_nafld_mtb/mtb_counts_g1pg2p_jrpca.txt",sep = "\t",row.names = FALSE,quote=FALSE)

masst_results<-fread("NAFLD_human_val/caussy_nafld_mtb/masst_human_to_mouse/comb_masst_results.txt")
mtb_sub<-mtb%>%filter(`row ID` %in% masst_results$human_feat)
write.table(mtb_sub,"NAFLD_human_val/caussy_nafld_mtb/mtb_counts_g1pg2p_jrpca_mousehits.txt",sep = "\t",row.names = FALSE,quote=FALSE)

comb_md<-m16S_md%>%
  dplyr::select(sample_name,groups,adv_fibrosis,traintest)%>%
  filter(!(sample_name %in% c("CIR22.001","FS.CIR2.001")))
write.table(comb_md,"NAFLD_human_val/caussy_nafld_16s/metadata_jrpca_g2p.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#plot correlation

m16s_key<-read_qza("NAFLD_human_val/caussy_nafld_16s/04.sklearn.gg2.asv.taxonomy.qza")$data%>%
  dplyr::rename(FeatureID=Feature.ID)%>%
  dplyr::mutate(name=paste(FeatureID,Taxon,sep=" "))%>%
  dplyr::select(FeatureID,name)

mz<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/quantification_table/quantification_table-00000.csv")%>%
  dplyr::select(1:2)%>%
  dplyr::rename(FeatureID=`row ID`,mz=`row m/z`)%>%
  dplyr::filter(FeatureID %in% masst_results$human_feat)

mtb_key<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/DB_analogresult/e902b4b271a04788b8bc9a1f25b0afd3.tsv")%>%
  dplyr::rename(FeatureID=`#Scan#`)%>%
  right_join(.,mz,by="FeatureID")%>%
  dplyr::mutate(name=ifelse(is.na(Compound_Name),paste("mtb",FeatureID,"m/z",signif(mz,5),sep=" "),paste("mtb",FeatureID,"m/z",signif(mz,5),Compound_Name,sep=" ")),
                FeatureID=as.character(FeatureID))%>%
  filter(!is.na(Compound_Name))%>%
  dplyr::select(FeatureID,name)

mtb_sub<-mtb%>%filter(`row ID` %in% mtb_key$FeatureID)
write.table(mtb_sub,"NAFLD_human_val/caussy_nafld_mtb/mtb_counts_g1pg2p_jrpca_annot.txt",sep = "\t",row.names = FALSE,quote=FALSE)


key<-rbind(m16s_key,mtb_key,fill=TRUE)

corr_n<-fread("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_lachnohits/correlation_table/Correlation.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=featureid)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,key,by="FeatureID")%>%
  filter(!is.na(name))%>%
  dplyr::select(c(3513,3508:3512))%>%
  column_to_rownames("name")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,key,by="FeatureID")%>%
  filter(!is.na(name))%>%
  dplyr::select(c(30,2:24))%>%
  column_to_rownames("name")

write.table(corr_n, "/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_lachnohits/correlation_table/Correlation_clean.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(corr_n, "/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_lachnohits/correlation_table/Correlation_clean_annot.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

corr_n<-fread("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_lachnohits_annot/correlation_table/Correlation.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=featureid)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,key,by="FeatureID")%>%
  filter(!is.na(name))%>%
  dplyr::select(c(30,25:29))%>%
  column_to_rownames("name")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,key,by="FeatureID")%>%
  filter(!is.na(name))%>%
  dplyr::select(c(30,2:24))%>%
  column_to_rownames("name")

#write.table(corr_n, "/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_lachnohits/correlation_table/Correlation_clean.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(corr_n, "/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_lachnohits_annot//correlation_table/Correlation_clean_annot.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

corr_n<-fread("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_lachnohits_mouse/correlation_table/Correlation.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=featureid)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,key,by="FeatureID")%>%
  filter(!is.na(name))%>%
  dplyr::select(c(83,78:82))%>%
  column_to_rownames("name")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,key,by="FeatureID")%>%
  filter(!is.na(name))%>%
  dplyr::select(c(30,2:24))%>%
  column_to_rownames("name")

write.table(corr_n, "/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_lachnohits_mouse/correlation_table/Correlation_clean_annot.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

p<-pheatmap(corr_n,color= inferno(10))
ggsave("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_lachnohits/correlation_table/SFR24_0323_heatmap_annot.pdf", p,width =30 , height = 26)
ggsave("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_lachnohits_annot/correlation_table/SFR24_0323_heatmap_annot.pdf", p,width = 30, height = 26)
ggsave("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_lachnohits_mouse/correlation_table/SFR24_0612_heatmap_annot.pdf", p,width = 30, height = 26)

p<-pheatmap(corr_n, annotation_row = m16s_key_ext,annotation_col = mtb_key_ext,
            annotation_colors = ann_colors,color= inferno(10), main="Wk12")



