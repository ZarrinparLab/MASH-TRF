setwd("/mnt/zarrinpar/scratch/sfloresr/NASH_KF")

library(tidyverse)
library(data.table)
library(ggpubr)
library(ggvenn)
library(gplots)
library(viridis)
library(RColorBrewer)
library(BioVenn) 
library(VennDiagram)
library("qiime2R")
library(ggalluvial)
library(ggh4x)
library(pheatmap)

#16S validation on human NAFLD dataset
##################################################################

taxa<-read_qza("NAFLD_human_val/caussy_nafld_16s/04.sklearn.gg2.asv.taxonomy.qza")$data%>%
  dplyr::rename(FeatureID=Feature.ID)
write.table(taxa,"NAFLD_human_val/taxonomy.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("NAFLD_human_val/caussy_nafld_16s//NAFLD_metadata.tsv")%>%
  filter(groups=="g1p"|groups=="g3p")%>%
  mutate(sample_name = sub(".*?\\.", "", sample_name))%>%
  mutate(traintest=ifelse(sample_name %in% data_mod$sample_name, "test","train"))
write.table(md,"NAFLD_human_val/caussy_nafld_16s/NAFLD_metadata_clean.tsv",sep = "\t",row.names = FALSE,quote=FALSE)

set.seed(1234)
data_mod <- md %>% group_by(groups) %>% sample_frac(.3)

m16s<-read_qza("NAFLD_human_val/caussy_nafld_16s/04.taxonomy_filtered.asv.counts.qza")$data%>%
  as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  dplyr::select(FeatureID,all_of(md$sample_name))

write.table(m16s,"NAFLD_human_val/caussy_nafld_16s/04.taxonomy_filtered.asv.counts.g1pg3p.txt",sep = "\t",row.names = FALSE,quote=FALSE)

##################################################################

bdm<-fread("NAFLD_human_val/caussy_nafld_16s/birdman/04.taxonomy_filtered.asv.counts.g1pg3p.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(groups, Treatment('g1p'))[T.g3p]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(groups, Treatment('g1p'))[T.g3p]_hdi`=gsub("[(]|[)]","",`C(groups, Treatment('g1p'))[T.g3p]_hdi`))%>%
  separate(`C(groups, Treatment('g1p'))[T.g3p]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
         #cred=ifelse(min>0.5|max< -0.5,"credible","not_credible"))%>%
         #cred=ifelse(min>1|max< -1,"credible","not_credible"))%>%
  left_join(.,taxa, by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon, sep=" "))%>%
  filter(cred=="credible")%>%
  dplyr::select(FeatureID, ratio, min, max, Name)%>%
  arrange(ratio)%>%
  mutate(cat=ifelse(ratio<0,"Non-MAFLD","MAFLD"))%>%
  mutate(cat=factor(cat,levels=c("Non-MAFLD","MAFLD")))

#write.table(bdm,"NAFLD_human_val/caussy_nafld_16s/birdman/birdman_justcred_g1pg3p_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
#write.table(bdm,"NAFLD_human_val/caussy_nafld_16s/birdman/birdman_justcred_g1pg3p_results_cred0.5.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
write.table(bdm,"NAFLD_human_val/caussy_nafld_16s/birdman/birdman_justcred_g1pg3p_results_cred1.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

bdm$Name <- factor(bdm$Name,levels = bdm$Name)

p<-ggplot(bdm, aes(x =Name , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle("Non-MAFLD vs. MAFLD")+
  coord_flip()

#ggsave("NAFLD_human_val/caussy_nafld_16s/birdman/SFR24_0213_g1pg3p.pdf", plot=p,height=15, width=25)
#ggsave("NAFLD_human_val/caussy_nafld_16s/birdman/SFR24_0213_g1pg3p_cred0.5.pdf", plot=p,height=10, width=25)
ggsave("NAFLD_human_val/caussy_nafld_16s/birdman/SFR24_0213_g1pg3p_cred1.pdf", plot=p,height=6, width=25)

ggplot(bdm, aes(x =Name , y = ratio, ymin = min, ymax = max, fill=cat)) + 
  geom_bar(stat="identity", alpha=0.7)+
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(title="Non-MAFLD vs. MAFLD")+
  coord_flip()
ggsave("NAFLD_human_val/caussy_nafld_16s/birdman/SFR24_0213_g1pg3p_cred1_colored.pdf",height=6, width=26)


sub_bdm<-bdm%>%filter(ratio< -5| ratio>5)
sub_bdm$Name <- factor(sub_bdm$Name,levels = sub_bdm$Name)

p<-ggplot(sub_bdm, aes(x =Name , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ggtitle("Non-MAFLD vs. MAFLD")+
  coord_flip()

ggsave("NAFLD_human_val/caussy_nafld_16s/birdman/SFR24_0213_g1pg3p_fold5.pdf", plot=p,height=8, width=25)

##################################################################
#is there a hit for the lachno or oscillo we say from our gut microbes here?

blast_oscillo<-fread("NAFLD_human_val/caussy_nafld_16s//WUMWNPBE114-Alignment_oscillo.txt")
blast_lachno<-fread("NAFLD_human_val/caussy_nafld_16s//WUN4T7BF114-Alignment_lachno.txt")

bdm_oscillo<-bdm%>%
  filter(FeatureID %in% blast_oscillo$`subject acc.ver`)%>%
  mutate(mash_hit="oscillo")

bdm_lachno<-bdm%>%
  filter(FeatureID %in% blast_lachno$`subject acc.ver`)%>%
  mutate(mash_hit="lachno")

combhits<-rbind(bdm_oscillo,bdm_lachno)

ggplot(combhits, aes(x =Name , y = ratio, ymin = min, ymax = max, fill=cat)) + 
  geom_bar(stat="identity", alpha=0.7)+
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(title="Non-MAFLD vs. MAFLD")+
  coord_flip()
ggsave("NAFLD_human_val/caussy_nafld_16s/birdman/SFR24_0214_g1pg3p_cred0_matchtomouselachnooscillo_colored.pdf",height=3.5, width=25.5)

#BDM params
# beta_prior: float = 5.0
# inv_disp_sd: float = 0.5

###################################################################
#create input file for joint-rpca

#16S
##metadata
m16S_md<-fread("NAFLD_human_val/caussy_nafld_16s/NAFLD_metadata_clean.tsv")

##quant_tab
m16S<-fread("NAFLD_human_val/caussy_nafld_16s/04.taxonomy_filtered.asv.counts.g1pg3p.txt")%>%
  dplyr::select(-CIR22.001,-FS.CIR2.001)
write.table(m16S,"NAFLD_human_val/caussy_nafld_16s/m16S_counts_g1pg3p_jrpca.txt",sep = "\t",row.names = FALSE,quote=FALSE)

m16s_sub<-m16S%>%filter(FeatureID %in% combhits$FeatureID)
write.table(m16s_sub,"NAFLD_human_val/caussy_nafld_16s/m16S_counts_g1pg3p_jrpca_microhits.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#mtb
##metadata
mtb_md<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/metadata/clean_fecal_metadata.tsv")%>%
   mutate(sampleid=gsub("\\-","\\.",sampleid))%>%
  filter(sampleid %in% m16S_md$sample_name)%>%
  filter(!(filename %in% c("TWDS002_RB6_01_29701.mzXML","TWCP-001_RC11_01_29719.mzXML")))

##quant_tab
mtb<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/quantification_table/quantification_table-G1PG3P.txt")%>%
  dplyr::select(-`TWDS002_RB6_01_29701.mzXML Peak area`,-`TWCP-001_RC11_01_29719.mzXML Peak area`)
colnames(mtb)<-c("FeatureID",mtb_md$sampleid)
write.table(mtb,"NAFLD_human_val/caussy_nafld_mtb/mtb_counts_g1pg3p_jrpca.txt",sep = "\t",row.names = FALSE,quote=FALSE)

mtb_sub<-mtb%>%filter(FeatureID %in% mtb_key$`#Scan#`)
write.table(mtb_sub,"NAFLD_human_val/caussy_nafld_mtb/mtb_counts_g1pg3p_jrpca_annot.txt",sep = "\t",row.names = FALSE,quote=FALSE)

masst_results<-fread("NAFLD_human_val/caussy_nafld_mtb/masst_human_to_mouse/comb_masst_results.txt")
mtb_sub<-mtb%>%filter(FeatureID %in% masst_results$human_feat)
write.table(mtb_sub,"NAFLD_human_val/caussy_nafld_mtb/mtb_counts_g1pg3p_jrpca_mousehits.txt",sep = "\t",row.names = FALSE,quote=FALSE)


comb_md<-m16S_md%>%
  dplyr::select(sample_name,groups,adv_fibrosis,traintest)%>%
  filter(!(sample_name %in% c("CIR22.001","FS.CIR2.001")))
write.table(comb_md,"NAFLD_human_val/caussy_nafld_16s/metadata_jrpca.txt",sep = "\t",row.names = FALSE,quote=FALSE)

# setdiff(colnames(mtb),colnames(m16S))
# intersect(colnames(m16S),colnames(mtb))

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

key<-rbind(m16s_key,mtb_key,fill=TRUE)

corr_n<-fread("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_microhits/correlation_table/Correlation.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=featureid)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  #left_join(.,m16s_key,by="FeatureID")%>%
  dplyr::select(c(3269,3254:3267))%>%
  column_to_rownames("name")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  #left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(1:3253))%>%
  column_to_rownames("FeatureID")

write.table(corr_n, "/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_microhits/correlation_table/Correlation_clean.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

corr_n<-fread("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_microhits_mouse/correlation_table/Correlation.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=featureid)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,key,by="FeatureID")%>%
  filter(!is.na(name))%>%
  dplyr::select(c(93,79:92))%>%
  #dplyr::select(c(548,534:547))%>%
  column_to_rownames("name")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,key,by="FeatureID")%>%
  filter(!is.na(name))%>%
  dplyr::select(c(39,2:24))%>%
  #dplyr::select(c(93,2:78))%>%
  column_to_rownames("name")

write.table(corr_n, "/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_microhits_annot/correlation_table/Correlation_clean.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(corr_n, "/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_microhits_annot/correlation_table/Correlation_clean_annot.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(corr_n, "/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_microhits_mouse/correlation_table/Correlation_clean_annot.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(corr_n, "/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_microhits_mouse/correlation_table/Correlation_clean_justannot.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

p<-pheatmap(corr_n,color= inferno(10))
#ggsave("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_microhits/correlation_table/SFR24_0323_heatmap_annot.pdf", p,width = 30, height = 5)
ggsave("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_microhits_annot/correlation_table/SFR24_0323_heatmap_annot.pdf", p,width = 49, height = 40)
ggsave("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_microhits_mouse/correlation_table/SFR24_0612_heatmap_annot.pdf", p,width = 46, height = 29)
ggsave("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/NAFLD_human_val/caussy_nafld_16s/result_microhits_mouse/correlation_table/SFR24_0612_heatmap_justannot.pdf", p,width = 30, height = 27)

p<-pheatmap(corr_n, annotation_row = m16s_key_ext,annotation_col = mtb_key_ext,
            annotation_colors = ann_colors,color= inferno(10), main="Wk12")



