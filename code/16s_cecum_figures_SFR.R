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

#16S Reanalysis of STAM-HCC data (cecum 16s)
##################################################################
#clean metadata

mdx<-fread("16s/ileum/metadata_cln_ileum_bdm.txt")%>%
  dplyr::select(host_subject_id, NASH_category)

md<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/cecum/data/raw-data-pre-10956-artifact-119356/metadata_study13785_from_prep_10956.txt")%>%
  filter(sample_type=="cecum tissue"& host_age==12)%>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  left_join(.,mdx,by="host_subject_id")%>%
  dplyr::rename(sample_name=`sample id`)%>%
  mutate(condition_ZT=paste(condition,collection_timepoint, sep="_"),
         NASH_category=ifelse(condition=="NA", "not_applicable", NASH_category))%>%
  mutate(NASH_category=ifelse(is.na(NASH_category), "unknown", NASH_category))

write.table(md,"16s/cecum/metadata_cln.txt",sep = "\t",row.names = FALSE, quote=FALSE)  

grade <- fread("files_from_KF/STAM_TRF/Histology/nas_scoring_analysis_nash_jingjing.csv", header=TRUE)%>%
  dplyr::select(mouse_id,steatosis_grade,fibrosis_stage)%>%
  dplyr::rename(host_subject_id=mouse_id)%>%
  mutate(steatosis_grade_new=ifelse(steatosis_grade>0, "1+","0"),
         fibrosis_stage_new=ifelse(fibrosis_stage>0, "1+","0"))

md<-md%>%left_join(.,grade, by="host_subject_id")%>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  mutate(steatosis_grade_new=ifelse(condition=="NA","none",steatosis_grade_new),
         fibrosis_stage_new=ifelse(condition=="NA","none",fibrosis_stage_new))

write.table(md,"16s/cecum/metadata_cln_addmoreNASHcat.txt",sep = "\t",row.names = FALSE, quote=FALSE)
##################################################################

#plot RPCA of all data

ord <- read_qza("16s/cecum/rpca_results/rpca_results_NASH_all/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/cecum/rpca_results/rpca_results_NASH_all/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/cecum/rpca_results/rpca_results_NASH_all/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/cecum/metadata_cln.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(collection_timepoint=="ZT1","light","dark"))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=phase)) +
  geom_point(alpha=1.0, size=2.5) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = condition_ZT))+
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH cecum 16S (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/cecum/rpca_results/rpca_results_NASH_all/SFR23_1012_NASH_cecum16s_RPCA.pdf", plot=p,height=4, width=4)

#plot ZT1 RPCA

ord <- read_qza("16s/cecum/rpca_results/rpca_results_NASH_ZT1/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/cecum/rpca_results/rpca_results_NASH_ZT1//sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/cecum/rpca_results/rpca_results_NASH_ZT1/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/cecum/metadata_cln.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(collection_timepoint=="ZT1","light","dark"))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=phase)) +
  geom_point(alpha=1.0, size=2.5, stroke=1.5) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = phase))+
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH cecum 16S ZT1 (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/cecum/rpca_results/rpca_results_NASH_ZT1/SFR23_1012_NASH_cecum16s_ZT1_RPCA.pdf", plot=p,height=4, width=4)

#plot ZT13 RPCA

ord <- read_qza("16s/cecum/rpca_results/rpca_results_NASH_ZT13/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/cecum/rpca_results/rpca_results_NASH_ZT13/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/cecum/rpca_results/rpca_results_NASH_ZT13/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/cecum/metadata_cln.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(collection_timepoint=="ZT1","light","dark"))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=phase)) +
  geom_point(alpha=1.0, size=2.5, stroke=1.5) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = phase))+
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_shape_manual(values=c(16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH cecum 16S ZT13 (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/cecum/rpca_results/rpca_results_NASH_ZT13/SFR23_1012_NASH_cecum16s_ZT13_RPCA.pdf", plot=p,height=4, width=4)


