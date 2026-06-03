setwd("~/Notebooks/sfloresr/MASH-TRF/MASHomics/")

library(tidyverse)
library(data.table)
library(ggpubr)
library("qiime2R")
library(cowplot)

##################################################################
#inputs
mash_metadata<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/NASH_stool_metabolomics_metadata_new.txt"
results_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/stool_mtb/"

#########################################################################
#functions

rpca_dat<-function(ord){
  
  rpca<-ord$data$Vectors %>%
    dplyr::select(SampleID, PC1, PC2, PC3)%>%
    left_join(md,by="SampleID")%>%
    mutate(ATTRIBUTE_disease_stage=factor(ATTRIBUTE_disease_stage,levels=c("TRF","NASH")),
           NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH")))
  return(rpca)
}

rpca_plot_disease_timediff<-function(rpca, disease, time, color_num) {
  if (time == "collection_disease_stage") {
    if (disease == "NASH_category") {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = NASH_category, shape = ATTRIBUTE_disease_stage))
    } else if (disease == "fibrosis_stage") {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = fibrosis_stage, shape = ATTRIBUTE_disease_stage))
    } else {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = steatosis_grade, shape = ATTRIBUTE_disease_stage))
    }
  } else {
    if (disease == "NASH_category") {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = NASH_category, shape = ATTRIBUTE_timepoint))
    } else if (disease == "fibrosis_stage") {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = fibrosis_stage, shape = ATTRIBUTE_timepoint))
    } else {
      p <- rpca %>% ggplot(aes(x = PC1, y = PC2, color = steatosis_grade, shape = ATTRIBUTE_timepoint))
    }
  }
  
  p <- p + geom_point(alpha = 1.0) +
    theme_classic() + theme(legend.position = "top") +
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

#########################################################################
#RPCA by condition and collection time--Fig. 3A
##Need to run 02-run_rpca_stool16mtb.sh first

ord <- read_qza(paste0(results_dir,"smtb_rpca_results_NASH_all/ordination.qza"))

samp_ord<-ord$data$Vectors
write.table(samp_ord,paste0(results_dir,"smtb_rpca_results_NASH_all/sample_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,paste0(results_dir,"smtb_rpca_results_NASH_all/feature_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread(mash_metadata)%>%
  mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition),
         age_in_weeks=ifelse(ATTRIBUTE_disease_stage=="TRF", "Wk8","Wk12"),
         stage_condition=paste(age_in_weeks,ATTRIBUTE_condition,sep="_"),
         mashZT=paste(ATTRIBUTE_timepoint, NASH_category,sep="_"))%>%
  dplyr::rename(SampleID=sample_name)


rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  left_join(md,by="SampleID")%>%
  mutate(phase=ifelse(ATTRIBUTE_timepoint=="ZT1","light","dark"))%>%
  mutate(ATTRIBUTE_condition=factor(ATTRIBUTE_condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         stage_condition=factor(stage_condition,levels=c("Wk8_NA","Wk12_NA","Wk8_FA","Wk12_FA","Wk8_FT","Wk12_FT")),
         ATTRIBUTE_disease_stage=factor(ATTRIBUTE_disease_stage,levels=c("TRF","NASH")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=stage_condition, shape=phase)) +
  geom_point(alpha=1.0, size=2.5) + 
  theme_classic() +
  scale_color_manual(values=c("#56B4E9","#0072B2","#E69F00","#D55E00","limegreen","darkgreen"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH stool mtb")+ 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave(paste0(results_dir,"smtb_rpca_results_NASH_all/SFR24_0125_NASH_stoolmtb_RPCA_sepwk.pdf"), plot=p,height=4, width=4)

##add box plot to axes

cond_rpca <- ggplot(rpca, aes(x =ATTRIBUTE_condition, y = PC2,fill=ATTRIBUTE_condition)) +
  geom_boxplot(alpha=0.3) + geom_point(colour="black", aes(fill=ATTRIBUTE_condition),
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

pairwise.t.test(rpca$PC2, rpca$ATTRIBUTE_condition, p.adjust.method = "fdr")

diet_rpca <- ggplot(rpca, aes(x =ATTRIBUTE_diet, y = PC1))+
  geom_boxplot(alpha=0.3) + geom_point(fill="black",colour="white",shape=21,position=position_jitter(width=0.1,height = 0.1),size=3) +
  scale_x_discrete(limits=rev)+
  scale_y_continuous(breaks=c(-0.1,0.,0.1))+
  labs(x = NA, y = "PC1") +
  # theme_void()
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+ coord_flip()

pairwise.t.test(rpca$PC1, rpca$ATTRIBUTE_diet, p.adjust.method = "fdr")
# NCD   
# HFD <2e-16

plot_a<-plot_grid(NULL,cond_rpca, nrow=3,rel_heights = c(2, 25,0.5))
plot_b<-plot_grid(plot_a,p, rel_widths = c(1, 3))
plot_c<-plot_grid(NULL,diet_rpca, rel_widths = c(2.8,8))
final_plot <- plot_grid(plot_b,plot_c,nrow=2,rel_heights  = c(4, 1))

ggsave(paste0(results_dir,"smtb_rpca_results_NASH_all/SFR24_0125_NASH_stoolmtb_RPCA_sepwk_addbxplt.pdf"), plot=final_plot,height=4.7, width=5.7)

##################################################################
#RPCA by condition and collection time for FA and FT only--Fig. 3B
##Need to run 02-run_rpca_stool16S.sh first

##Wk 4 (overall week 8), left

ord <- read_qza(paste0(results_dir,"smtb_rpca_results_NASH_Wk8_FAFT/ordination.qza"))

samp_ord<-ord$data$Vectors
write.table(samp_ord,paste0(results_dir,"smtb_rpca_results_NASH_Wk8_FAFT/sample_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,paste0(results_dir,"smtb_rpca_results_NASH_Wk8_FAFT/feature_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  left_join(md,by="SampleID")%>%
  mutate(phase=ifelse(ATTRIBUTE_timepoint=="ZT1","light","dark"))%>%
  mutate(ATTRIBUTE_condition=factor(ATTRIBUTE_condition,levels=c("FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","NA")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_condition, shape=ZT)) +
  geom_point(alpha=1.0, size=4) + 
  theme_classic() +
  scale_color_manual(values=c("#D55E00","#009E73"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH stool mtb Wk8")+ 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave(paste0(results_dir,"smtb_rpca_results_NASH_Wk8_FAFT/SFR24_0125_NASH_stoolmtb_RPCA_sepLD.pdf"), plot=p,height=4, width=4)

##Wk 7 (overall week 12), right
ord <- read_qza(paste0(results_dir,"smtb_rpca_results_NASH_Wk12_FAFT/ordination.qza"))

samp_ord<-ord$data$Vectors
write.table(samp_ord,paste0(results_dir,"smtb_rpca_results_NASH_Wk12_FAFT/sample_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,paste0(results_dir,"smtb_rpca_results_NASH_Wk12_FAFT/feature_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  left_join(md,by="SampleID")%>%
  mutate(phase=ifelse(ATTRIBUTE_timepoint=="ZT1","light","dark"))%>%
  mutate(ATTRIBUTE_condition=factor(ATTRIBUTE_condition,levels=c("FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","NA")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=ATTRIBUTE_condition, shape=ZT)) +
  geom_point(alpha=1.0, size=4) + 
  theme_classic() +
  scale_color_manual(values=c("#D55E00","#009E73"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH stool mtb Wk12")+ 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave(paste0(results_dir,"smtb_rpca_results_NASH_Wk12_FAFT/SFR24_0125_NASH_stoolmtb_RPCA_sepLD.pdf"), plot=p,height=4, width=4)
