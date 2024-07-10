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

#Transcriptomics Reanalysis of STAM-HCC data (cecum)
##################################################################

nash_score<-fread("transcriptomics/md_with_nash_score.tsv")%>%
  dplyr::select(sample_id,NASH_category)
md_cec<-fread("transcriptomics/metadata.stam.noblanks.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  filter(sample_type=="cecum" & disease_state=="nash")%>%
  left_join(.,nash_score,by="sample_id")%>%
  dplyr::rename(sampleid=sample_id)%>%
  mutate(NASH_category=ifelse(condition=="NA","not_applicable",NASH_category))

write.table(md_cec,"transcriptomics/cecum_md_with_nash_score.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

grade <- fread("files_from_KF/STAM_TRF/Histology/nas_scoring_analysis_nash_jingjing.csv", header=TRUE)%>%
  dplyr::select(mouse_id,steatosis_grade,fibrosis_stage)%>%
  dplyr::rename(host_subject_id=mouse_id)%>%
  mutate(steatosis_grade_new=ifelse(steatosis_grade>0, "1+","0"),
         fibrosis_stage_new=ifelse(fibrosis_stage>0, "1+","0"))

md_cec<-md_cec%>%left_join(.,grade, by="host_subject_id")%>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  mutate(steatosis_grade_new=ifelse(condition=="NA","none",steatosis_grade_new),
         fibrosis_stage_new=ifelse(condition=="NA","none",fibrosis_stage_new))

write.table(md_cec,"transcriptomics/cecum_md_with_nash_score_fibste.tsv",sep = "\t",row.names = FALSE, quote=FALSE)


md_cec_ZT1<-md_cec%>%filter(timepoint=="ZT1")
md_cec_ZT13<-md_cec%>%filter(timepoint=="ZT13")
md_cec_NA<-md_cec%>%filter(condition=="NA")
md_cec_FA<-md_cec%>%filter(condition=="FA")
md_cec_FT<-md_cec%>%filter(condition=="FT")

#convert data file
dat<-fread("transcriptomics/data_files/cecum.nash.stringtie.gene_count_matrix.csv")
write.table(dat,"transcriptomics/data_files/cecum.nash.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_ZT1<-dat%>%dplyr::select(gene_id,md_cec_ZT1$sampleid)
write.table(dat_ZT1,"transcriptomics/data_files/cecum.nash.ZT1.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_ZT13<-dat%>%dplyr::select(gene_id,md_cec_ZT13$sampleid)
write.table(dat_ZT13,"transcriptomics/data_files/cecum.nash.ZT13.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_NA<-dat%>%dplyr::select(gene_id,md_cec_NA$sampleid)
write.table(dat_NA,"transcriptomics/data_files/cecum.nash.NA.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_FA<-dat%>%dplyr::select(gene_id,md_cec_FA$sampleid)
write.table(dat_FA,"transcriptomics/data_files/cecum.nash.FA.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_FT<-dat%>%dplyr::select(gene_id,md_cec_FT$sampleid)
write.table(dat_FT,"transcriptomics/data_files/cecum.nash.FT.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

##################################################################

ord <- read_qza("transcriptomics/rpca_results_cecum/rpca_results_NASH_all/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"transcriptomics/rpca_results_cecum/rpca_results_NASH_all/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"transcriptomics/rpca_results_cecum/rpca_results_NASH_all/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("transcriptomics/cecum_md_with_nash_score_fibste.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(sample_name=sampleid)

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(timepoint=="ZT1","light","dark"))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","not_applicable")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=phase)) +
  geom_point(alpha=1.0, size=2.5) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = condition_ZT))+
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  #scale_shape_manual(values=c(16,17,3)) +
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH cecum RNA (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("transcriptomics/rpca_results_cecum/rpca_results_NASH_all/SFR23_1024_NASH_cecumRNA_RPCA.pdf", plot=p,height=4, width=4)
ggsave("transcriptomics/rpca_results_cecum/rpca_results_NASH_all/SFR23_1024_NASH_cecumRNA_NASHcat_RPCA.pdf", plot=p,height=4, width=4)
