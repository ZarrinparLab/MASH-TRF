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

#16S Reanalysis of STAM-HCC data (ileum 16s)

##################################################################
#clean metadata to include a category of Pre-NASH (8 weeks) and NASH (12 weeks)

md<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/ileum/data/raw-data-pre-12728-artifact-143225/metadata_ileum_study13785_from_prep_12728.txt")%>%
  filter(empo_1!="Control" & host_age==12)%>% 
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  mutate(condition_ZT=paste(condition,collection_timepoint, sep="_"),
         NASH_category=ifelse(condition=="NA", "not_applicable", NASH_category))

write.table(md,"16s/ileum/metadata_cln.txt",sep = "\t",row.names = FALSE, quote=FALSE)  

grade <- fread("files_from_KF/STAM_TRF/Histology/nas_scoring_analysis_nash_jingjing.csv", header=TRUE)%>%
  dplyr::select(mouse_id,steatosis_grade,fibrosis_stage)%>%
  dplyr::rename(host_subject_id=mouse_id)%>%
  mutate(steatosis_grade_new=ifelse(steatosis_grade>0, "1+","0"),
         fibrosis_stage_new=ifelse(fibrosis_stage>0, "1+","0"))

md<-md%>%left_join(.,grade, by="host_subject_id")%>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  mutate(steatosis_grade_new=ifelse(condition=="NA","none",steatosis_grade_new),
         fibrosis_stage_new=ifelse(condition=="NA","none",fibrosis_stage_new))

write.table(md,"16s/ileum/metadata_cln_addmoreNASHcat.txt",sep = "\t",row.names = FALSE, quote=FALSE)

md_bdm<-md%>%dplyr::rename(sample_name=sampleid)
write.table(md_bdm,"16s/ileum/metadata_cln_ileum_bdm.txt",sep = "\t",row.names = FALSE, quote=FALSE)  
##################################################################
#plot RPCA of all data

ord <- read_qza("16s/ileum/rpca_results/rpca_results_NASH_all/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/ileum/rpca_results/rpca_results_NASH_all/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/ileum/rpca_results/rpca_results_NASH_all/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/ileum/metadata_cln_ileum_bdm.txt")%>%
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
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH ileum 16S (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/ileum/rpca_results/rpca_results_NASH_all/SFR23_1012_NASH_ileum16s_RPCA.pdf", plot=p,height=4, width=4)

#plot NA RPCA

ord <- read_qza("16s/ileum/rpca_results/rpca_results_NASH_NA/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/ileum/rpca_results/rpca_results_NASH_NA/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/ileum/rpca_results/rpca_results_NASH_NA/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/ileum/metadata_cln_ileum_bdm.txt")%>%
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
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = phase))+
  scale_color_manual(values=c("#0072B2"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH ileum 16S NA (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/ileum/rpca_results/rpca_results_NASH_NA/SFR23_1012_NASH_ileum16s_NA_RPCA.pdf", plot=p,height=4, width=4)

#plot FA RPCA

ord <- read_qza("16s/ileum/rpca_results/rpca_results_NASH_FA/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/ileum/rpca_results/rpca_results_NASH_FA/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/ileum/rpca_results/rpca_results_NASH_FA/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/ileum/metadata_cln_ileum_bdm.txt")%>%
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
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = phase))+
  scale_color_manual(values=c("#D55E00"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH ileum 16S FA (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/ileum/rpca_results/rpca_results_NASH_FA/SFR23_1012_NASH_ileum16s_FA_RPCA.pdf", plot=p,height=4, width=4)

#plot FT RPCA

ord <- read_qza("16s/ileum/rpca_results/rpca_results_NASH_FT/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/ileum/rpca_results/rpca_results_NASH_FT/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/ileum/rpca_results/rpca_results_NASH_FT/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/ileum/metadata_cln_ileum_bdm.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(collection_timepoint=="ZT1","light","dark"))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         NASH_category=factor(NASH_category, levels=c("Non_NASH","NASH")))

p<-rpca %>%
  #ggplot(aes(x=PC1, y=PC2, color=condition, shape=phase)) +
  ggplot(aes(x=PC1, y=PC2, color=NASH_category, shape=phase)) +
  geom_point(alpha=1.0, size=2.5) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = phase))+
  #scale_color_manual(values=c("#009E73"))+
  scale_color_manual(values=c("#0000a7","#c1272d"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
  #labs(color="NASH Category",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+
  ggtitle("NASH ileum 16S FT (Wk 12)")+ 
  theme(plot.title = element_text(face = "bold"))
#ggsave("16s/ileum/rpca_results/rpca_results_NASH_FT/SFR23_1012_NASH_ileum16s_FT_RPCA.pdf", plot=p,height=4, width=4)
ggsave("16s/ileum/rpca_results/rpca_results_NASH_FT/SFR23_1012_NASH_ileum16s_FT_nonnashnash_RPCA.pdf", plot=p,height=4, width=4)
##################################################################

#plot BIRDMAn results

##non-NASH VS. NASH ZT1 FA (RPCA showed no diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

ZT1FA<-fread("16s/stool/birdman_results/NASH_ZT1_FA_taxonomy_filtered.asv.counts.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`=gsub("[(]|[)]","",`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`))%>%
  separate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  left_join(.,annot, by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon, sep=" "))%>%
  #filter(FeatureID %in% ZT13FA$FeatureID)%>%
  filter(cred=="credible")%>%
  dplyr::select(FeatureID, ratio, min, max, Name)%>%
  arrange(ratio)%>%
  mutate(ZT_time="ZT1")

write.table(ZT1FA,"16s/stool/birdman_results/birdman_justcred_ZT1FA_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

ZT1FA$Name <- factor(ZT1FA$Name,levels = ZT1FA$Name)

ggplot(ZT1FA, aes(x =Name , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title="FA ZT1\nnon-NASH vs. NASH")+
  coord_flip()

ggsave("16s/stool/birdman_results/birdman_justcred_ZT1FA_nonNASHNASH.pdf",height=6, width=13.5)

##non-NASH VS. NASH ZT13 FA (RPCA showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

ZT13FA<-fread("16s/stool/birdman_results/NASH_ZT13_FA_taxonomy_filtered.asv.counts.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`=gsub("[(]|[)]","",`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`))%>%
  separate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  left_join(.,annot, by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon, sep=" "))%>%
  #filter(cred=="credible")%>%
  filter(FeatureID %in% ZT1FA$FeatureID)%>%
  dplyr::select(FeatureID, ratio, min, max, Name)%>%
  arrange(ratio)%>%
  mutate(ZT_time="ZT13")

write.table(ZT13FA,"16s/stool/birdman_results/birdman_justcred_ZT13FA_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

ZT13FA$Name <- factor(ZT13FA$Name,levels = ZT13FA$Name)

ggplot(ZT13FA, aes(x =Name , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title="FA ZT13\nnon-NASH vs. NASH")+
  coord_flip()

ggsave("16s/stool/birdman_results/birdman_justcred_ZT13FA_nonNASHNASH.pdf",height=6, width=14.5)


###plot ZT1 and ZT13 on the same plot
NASH_FA<-rbind(ZT1FA,ZT13FA)

#write.table(NASH_FA,"16s/stool/birdman_results/birdman_justcred_ZT1ZT13FA_ZT13ord_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
write.table(NASH_FA,"16s/stool/birdman_results/birdman_justcred_ZT1ZT13FA_ZT1ord_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

#NASH_FA$Name <- factor(NASH_FA$Name,levels = ZT13FA$Name)
NASH_FA$Name <- factor(NASH_FA$Name,levels = ZT1FA$Name)

ggplot(NASH_FA, aes(x =Name , y = ratio, ymin = min, ymax = max, color=ZT_time, group=ZT_time)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +scale_color_manual(values=c("gray70","gray10")) +
  labs(title="FA non-NASH vs. NASH")

#ggsave("16s/stool/birdman_results/birdman_justcred_FA_ordZT13_nonNASHNASH.pdf",height=6, width=14.5)
ggsave("16s/stool/birdman_results/birdman_justcred_FA_ordZT1_nonNASHNASH.pdf",height=6, width=14.5)

##non-NASH VS. NASH ZT1 FT (RPCA showed borderline diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

ZT1FT<-fread("16s/stool/birdman_results/NASH_ZT1_FT_taxonomy_filtered.asv.counts.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`=gsub("[(]|[)]","",`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`))%>%
  separate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  left_join(.,annot, by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon, sep=" "))%>%
  filter(cred=="credible")%>%
  #filter(FeatureID %in% ZT13FT$FeatureID)%>%
  dplyr::select(FeatureID, ratio, min, max, Name)%>%
  arrange(ratio)%>%
  mutate(ZT_time="ZT1")

write.table(ZT1FT,"16s/stool/birdman_results/birdman_justcred_ZT1FT_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

ZT1FT$Name <- factor(ZT1FT$Name,levels = ZT1FT$Name)

ggplot(ZT1FT, aes(x =Name , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title="FT ZT1\nnon-NASH vs. NASH")+
  coord_flip()

ggsave("16s/stool/birdman_results/birdman_justcred_ZT1FT_nonNASHNASH.pdf",height=6, width=14.5)

##non-NASH VS. NASH ZT13 FT (RPCA showed no diff--8 vs 12 wk diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

ZT13FT<-fread("16s/stool/birdman_results/NASH_ZT13_FT_taxonomy_filtered.asv.counts.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`=gsub("[(]|[)]","",`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`))%>%
  separate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  left_join(.,annot, by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon, sep=" "))%>%
  filter(cred=="credible")%>%
  filter(FeatureID %in% ZT1FT$FeatureID)%>%
  dplyr::select(FeatureID, ratio, min, max, Name)%>%
  arrange(ratio)%>%
  mutate(ZT_time="ZT13")

write.table(ZT13FT,"16s/stool/birdman_results/birdman_justcred_ZT13FT_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

ZT13FT$Name <- factor(ZT13FT$Name,levels = ZT13FT$Name)

ggplot(ZT13FT, aes(x =Name , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title="FT ZT13\nnon-NASH vs. NASH")+
  coord_flip()

ggsave("16s/stool/birdman_results/birdman_justcred_ZT13FT_nonNASHNASH.pdf",height=4, width=14.5)

###plot ZT1 and ZT13 on the same plot
NASH_FT<-rbind(ZT1FT,ZT13FT)

#write.table(NASH_FT,"16s/stool/birdman_results/birdman_justcred_ZT1ZT13FT_ZT13ord_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
write.table(NASH_FT,"16s/stool/birdman_results/birdman_justcred_ZT1ZT13FT_ZT1ord_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

#NASH_FT$Name <- factor(NASH_FT$Name,levels = ZT13FT$Name)
NASH_FT$Name <- factor(NASH_FT$Name,levels = ZT1FT$Name)

ggplot(NASH_FT, aes(x =Name , y = ratio, ymin = min, ymax = max, color=ZT_time, group=ZT_time)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +scale_color_manual(values=c("gray70","gray10")) +
  labs(title="FT non-NASH vs. NASH")

#ggsave("16s/stool/birdman_results/birdman_justcred_FT_ordZT13_nonNASHNASH.pdf",height=5, width=14.5)
ggsave("16s/stool/birdman_results/birdman_justcred_FT_ordZT1_nonNASHNASH.pdf",height=6, width=14.5)



##8 week VS. 12 week ZT13 FT (RPCA showed no diff--8 vs 12 wk diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

ZT13FT<-fread("16s/stool/birdman_results/NASH_ZT13_FT_wkcompar_taxonomy_filtered.asv.counts.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(collection_disease_stage, Treatment('TRF'))[T.NASH]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(collection_disease_stage, Treatment('TRF'))[T.NASH]_hdi`=gsub("[(]|[)]","",`C(collection_disease_stage, Treatment('TRF'))[T.NASH]_hdi`))%>%
  separate(`C(collection_disease_stage, Treatment('TRF'))[T.NASH]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  left_join(.,annot, by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon, sep=" "))%>%
  filter(cred=="credible")%>%
  dplyr::select(FeatureID, ratio, min, max, Name)%>%
  arrange(ratio)

write.table(ZT13FT,"16s/stool/birdman_results/birdman_justcred_ZT13FT_8wkv12wk_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

ZT13FT$Name <- factor(ZT13FT$Name,levels = ZT13FT$Name)

ggplot(ZT13FT, aes(x =Name , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title="FT ZT13\n8week vs. 12 week")+
  coord_flip()

ggsave("16s/stool/birdman_results/birdman_justcred_ZT13FT_8wkv12wk.pdf",height=6, width=14.5)

##################################################################

#comparing overlaps 

##ZT1 vs ZT13 for FA

ZT1FA<-fread("16s/stool/birdman_results/birdman_justcred_ZT1FA_nonNASHNASH_results.tsv")
ZT13FA<-fread("16s/stool/birdman_results/birdman_justcred_ZT13FA_nonNASHNASH_results.tsv")
  
list_venn <- list(ZT1 = ZT1FA$Name,
                  ZT13 = ZT13FA$Name)

p<-venn.diagram(list_venn, fill = c("white", "black"),height = 10,
                main="NASH FA",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_FA_nonNASHNASH_results.pdf")
grid.draw(p)
dev.off()

ZT1FT<-fread("16s/stool/birdman_results/birdman_justcred_ZT1FT_nonNASHNASH_results.tsv")
ZT13FT<-fread("16s/stool/birdman_results/birdman_justcred_ZT13FT_nonNASHNASH_results.tsv")

list_venn <- list(ZT1 = ZT1FT$Name,
                  ZT13 = ZT13FT$Name)

p<-venn.diagram(list_venn, fill = c("white", "black"),height = 10,
                main="NASH FT",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_FT_nonNASHNASH_results.pdf")
grid.draw(p)
dev.off()


list_venn <- list(ZT1_FA = ZT1FA$Name,
                  ZT13_FA = ZT13FA$Name,
                  ZT1_FT = ZT1FT$Name,
                  ZT13_FT = ZT13FT$Name)

p<-venn.diagram(list_venn,height = 10,
                main="NASH",
                width = 10,alpha = c(0.5, 0.5, 0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_nonNASHNASH_results.pdf")
grid.draw(p)
dev.off()
