setwd("~/Notebooks/sfloresr/MASH-TRF/MASHomics/")

library(tidyverse)
library(data.table)
library(ggpubr)
library("qiime2R")
library(cowplot)

##################################################################
#inputs
mash_metadata<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_16s/mashstool16s_preprocessed_20211020_ID_13785_gg2/metadata.txt"
NASH_histology<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/nas_scoring_analysis_nash.csv"
data_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_16s/"
results_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/stool_16s/"
##################################################################
#functions

rpca_dat<-function(ord){
  
  rpca<-ord$data$Vectors %>%
    dplyr::select(SampleID, PC1, PC2, PC3)%>%
    dplyr::rename(sample_name=SampleID)%>%
    left_join(md,by="sample_name")%>%
    mutate(collection_disease_stage=factor(collection_disease_stage,levels=c("TRF","NASH")),
           NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH")))
  return(rpca)
}

rpca_plot_disease_timediff<-function(rpca, disease, time, color_num) {
  if (time == "collection_disease_stage") {
    if (disease == "NASH_category") {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = NASH_category, shape = collection_disease_stage))
    } else if (disease == "fibrosis_stage") {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = fibrosis_stage, shape = collection_disease_stage))
    } else {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = steatosis_grade, shape = collection_disease_stage))
    }
  } else {
    if (disease == "NASH_category") {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = NASH_category, shape = collection_timepoint))
    } else if (disease == "fibrosis_stage") {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = fibrosis_stage, shape = collection_timepoint))
    } else {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = steatosis_grade, shape = collection_timepoint))
    }
  }
  
  p <- p + geom_point(alpha = 1.0) +
    theme_pubr() +
    scale_shape_manual(values = c(3, 16)) +
    labs(
      color = disease,
      x = paste("PC1 (", round(ord$data$ProportionExplained$PC1 * 100, digits = 2), "%)", sep = ""),
      y = paste("PC2 (", round(ord$data$ProportionExplained$PC2 * 100, digits = 2), "%)", sep = "")
    )
  
  if (color_num == 2) {
    p <- p + scale_color_manual(values = c("#0000a7", "#c1272d"))
  } else {
    p <- p + scale_color_manual(values = c("#648FFF", "#FFB000", "#DC267F"))
  }
  
  return(p)
}

##################################################################
#clean metadata

##include a category of Pre-NASH (8 weeks) and NASH (12 weeks)
md<-fread(mash_metadata)%>%
  filter(cohort=="NASH")%>% 
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  mutate(condition_ZT=paste(condition,collection_timepoint, sep="_"),
         NASH_category=ifelse(NASH_category=="N/A", "not_applicable", NASH_category))

write.table(md,paste0(data_dir,"s16s_metadata_cln.txt"),sep = "\t",row.names = FALSE, quote=FALSE)  

md<-md%>%mutate(NASH_category_new=ifelse(host_age==8,"pre-NASH",NASH_category))
write.table(md,paste0(data_dir,"s16s_metadata_cln_addNASHcat.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

##add fibrosis stage and steatosis grade categories to md
grade <- fread(NASH_histology, header=TRUE)%>%
  dplyr::select(mouse_id,steatosis_grade,fibrosis_stage)%>%
  dplyr::rename(host_subject_id=mouse_id)%>%
  mutate(steatosis_grade_new=ifelse(steatosis_grade>0, "1+","0"),
         fibrosis_stage_new=ifelse(fibrosis_stage>0, "1+","0"))

md<-md%>%left_join(.,grade, by="host_subject_id")%>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  mutate(steatosis_grade_new=ifelse(condition=="NA","none",steatosis_grade_new),
         fibrosis_stage_new=ifelse(condition=="NA","none",fibrosis_stage_new))

write.table(md,paste0(data_dir,"s16s_metadata_cln_addmoreNASHcat.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

##################################################################
#RPCA by condition and collection time--Fig. 1A
##Need to run 01-run_rpca_stool16S.sh first

ord <- read_qza(paste0(results_dir,"s16s_rpca_results_NASH_all_gg2/ordination.qza"))

samp_ord<-ord$data$Vectors
write.table(samp_ord,paste0(results_dir,"s16s_rpca_results_NASH_all_gg2/sample_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,paste0(results_dir,"s16s_rpca_results_NASH_all_gg2/feature_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread(paste0(data_dir,"s16s_metadata_cln.txt"))%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         age_in_weeks=ifelse(collection_disease_stage=="TRF", "Wk8","Wk12"),
         stage_condition=paste(age_in_weeks,condition,sep="_"),
         diet=ifelse(condition=="NA","NCD","HFD"))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(collection_timepoint=="ZT1","light","dark"))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         diet=factor(diet,levels=c("NCD","HFD")),
         phase=factor(phase,levels=c("light","dark")),
         stage_condition=factor(stage_condition,levels=c("Wk8_NA","Wk12_NA","Wk8_FA","Wk12_FA","Wk8_FT","Wk12_FT")),
         collection_disease_stage=factor(collection_disease_stage,levels=c("TRF","NASH")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=stage_condition, shape=phase)) +
  geom_point(alpha=1.0, size=2.5) +
  theme_pubr() +
  scale_color_manual(values=c("#56B4E9","#0072B2","#E69F00","#D55E00","limegreen","darkgreen"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH stool 16S")+ 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")
ggsave(paste0(results_dir,"s16s_rpca_results_NASH_all_gg2/SFR24_0209_NASH_stool16s_RPCA_sepwk.pdf"), plot=p,height=4, width=4)

##add box plot to axes

cond_rpca <- ggplot(rpca, aes(x =condition, y = PC2,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_point(colour="black", aes(fill=condition),
                                       position=position_jitter(width=0.1,height = 0.1),
                                       size=3,pch=21) + 
  labs(x = NA, y = "PC2") +
  scale_y_continuous(breaks=c(-0.2,-0.1,0.,0.1,0.2,0.3))+
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  # theme_void()
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

pairwise.t.test(rpca$PC2, rpca$condition, p.adjust.method = "fdr")
# NA      FA     
# FA 0.00248 -      
#   FT 0.49580 0.00054


diet_rpca <- ggplot(rpca, aes(x =diet, y = PC1))+
  geom_boxplot(alpha=0.3) + geom_point(fill="black",colour="white",shape=21,position=position_jitter(width=0.1,height = 0.1),size=3) +
  scale_x_discrete(limits=rev)+
  scale_y_continuous(breaks=c(-0.1,0.,0.1))+
  labs(x = NA, y = "PC1") +
  # theme_void()
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+ coord_flip()

pairwise.t.test(rpca$PC1, rpca$diet, p.adjust.method = "fdr")
# NCD   
# HFD <2e-16

plot_a<-plot_grid(NULL,cond_rpca, nrow=3,rel_heights = c(2, 25,0.5))
plot_b<-plot_grid(plot_a,p, rel_widths = c(1, 3))
plot_c<-plot_grid(NULL,diet_rpca, rel_widths = c(2.8,8))
final_plot <- plot_grid(plot_b,plot_c,nrow=2,rel_heights  = c(4, 1))

ggsave(paste0(results_dir,"s16s_rpca_results_NASH_all_gg2/SFR24_0508_NASH_stool16s_RPCA_sepwk_addbxplt.pdf"), plot=final_plot,height=4.7, width=5.7)

##################################################################
#RPCA by condition and collection time for FA and FT only--Fig. 1B
##Need to run 01-run_rpca_stool16S.sh first

##Wk 4 (overall week 8), left

ord <- read_qza(paste0(results_dir,"s16_rpca_results_NASH_8wk_FAFT_gg2/ordination.qza"))

samp_ord<-ord$data$Vectors
write.table(samp_ord,paste0(results_dir,"s16s_rpca_results_NASH_8wk_FAFT_gg2/sample_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,paste0(results_dir,"s16s_rpca_results_NASH_8wk_FAFT_gg2/feature_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(collection_timepoint=="ZT1","light","dark"))%>%
  mutate(ATTRIBUTE_condition=factor(condition,levels=c("FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","NA")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=phase)) +
  geom_point(alpha=1.0, size=4) + 
  theme_classic() +
  scale_color_manual(values=c("#D55E00","#009E73"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH stool 16S Wk8")+ 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave(paste0(results_dir,"s16s_rpca_results_NASH_8wk_FAFT_gg2/SFR24_0209_NASH_stool16s_RPCA_sepLD.pdf"), plot=p,height=4, width=4)

##Wk 7 (overall week 12), right

ord <- read_qza(paste0(results_dir,"s16s_rpca_results_NASH_12wk_FAFT_gg2/ordination.qza"))

samp_ord<-ord$data$Vectors
write.table(samp_ord,paste0(results_dir,"s16s_rpca_results_NASH_12wk_FAFT_gg2/sample_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,paste0(results_dir,"s16s_rpca_results_NASH_12wk_FAFT_gg2/feature_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(collection_timepoint=="ZT1","light","dark"))%>%
  mutate(ATTRIBUTE_condition=factor(condition,levels=c("FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","NA")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=phase)) +
  geom_point(alpha=1.0, size=4) + 
  theme_classic() +
  scale_color_manual(values=c("#D55E00","#009E73"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH stool 16S Wk12")+ 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave(paste0(results_dir,"s16s_rpca_results_NASH_12wk_FAFT_gg2/SFR24_0209_NASH_stool16s_RPCA_sepLD.pdf"), plot=p,height=4, width=4)
