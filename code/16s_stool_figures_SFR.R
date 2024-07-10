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
library(cowplot)
#16S Reanalysis of STAM-HCC data (stool 16s)

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

rpca_plot_disease<-function(rpca,disease,color_num){
  if(disease=="NASH_category"){
    p<-rpca %>%
      ggplot(aes(x=PC1, y=PC2, color=NASH_category))
  }
  else if(disease=="fibrosis_stage"){
    p<-rpca %>%
      ggplot(aes(x=PC1, y=PC2, color=fibrosis_stage))
  }
  else{
    p<-rpca %>%
      ggplot(aes(x=PC1, y=PC2, color=steatosis_grade))
  }

  p<-p+ geom_point(alpha=1.0) +
    #stat_ellipse(type = "t", linetype = 2,aes(group = NASH_category))+
    theme_pubr() +
    labs(color=disease,
         x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
         y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))
  
  if (color_num==2) {
    p<-p + scale_color_manual(values=c("#0000a7","#c1272d"))
  }
  else
    {p<-p+scale_color_manual(values=c("#648FFF","#FFB000","#DC267F"))
  }
  return(p)
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
    #stat_ellipse(type = "t", linetype = 2,aes(group = NASH_category))+
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

birdman_dat<-function(dat,annot){
  df<-dat%>%
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
    dplyr::select(FeatureID, ratio, min, max, Name)%>%
    arrange(ratio)
  return(df)
}

birdman_plot<-function(dat){
  dat$Name <- factor(dat$Name,levels = dat$Name)
  
  p<-ggplot(dat, aes(x =Name , y = ratio, ymin = min, ymax = max)) + 
    geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+#theme_pubr()+
    geom_pointrange(position = position_dodge(width = 0.8)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip()
  return(p)
}

venn_plt<-function(df1,df2,title){
  df1<-fread("16s/stool/birdman_results/birdman_justcred_wk8FAZT1_nonNASHNASH_results.tsv")
  df2<-fread("16s/stool/birdman_results/birdman_justcred_wk8FAZT13_nonNASHNASH_results.tsv")
  
  list_venn <- list(ZT1 = df1$Name,
                    ZT13 = df2$Name)
  
  p<-venn.diagram(list_venn, fill = c("white", "black"),height = 10,
                  main=title,
                  width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)
  
  return(p)
}

paired_plots_ind<-function(dat,microid){
  count_dat<-dat%>%
    #filter(FeatureID==microid & condition=="FA" & 
    filter(FeatureID==microid & condition=="FT" & 
            NASH_category!="not_applicable")%>%
    arrange(collection_timepoint,collection_disease_stage)%>%
    arrange(host_subject_id)%>%as.data.table()
  
  p<-ggplot(count_dat, aes(x =NASH_category , y = counts)) +
    geom_boxplot(aes(fill = NASH_category), alpha = .5) +
    #geom_line(aes(group = host_subject_id)) + 
    ggtitle(stringr::str_wrap(unique(count_dat$Name), width=30))+
    #geom_point(size = 2) + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=2)+
    theme_bw()+ theme(legend.position = "top")+
    scale_fill_manual(values=c("#0000a7","#c1272d"))+
    facet_wrap(~ collection_disease_stage)
  return(p)
}
##################################################################
#clean metadata to include a category of Pre-NASH (8 weeks) and NASH (12 weeks)
md<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785/metadata.txt")%>%
  filter(cohort=="NASH")%>% 
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  mutate(condition_ZT=paste(condition,collection_timepoint, sep="_"),
         NASH_category=ifelse(NASH_category=="N/A", "not_applicable", NASH_category))

write.table(md,"16s/stool/metadata_cln.txt",sep = "\t",row.names = FALSE, quote=FALSE)  

md<-md%>%mutate(NASH_category_new=ifelse(host_age==8,"pre-NASH",NASH_category))
write.table(md,"16s/stool/metadata_cln_addNASHcat.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#add fibrosis stage and steatosis grade categories to md

grade <- fread("files_from_KF/STAM_TRF/Histology/nas_scoring_analysis_nash_jingjing.csv", header=TRUE)%>%
  dplyr::select(mouse_id,steatosis_grade,fibrosis_stage)%>%
  dplyr::rename(host_subject_id=mouse_id)%>%
  mutate(steatosis_grade_new=ifelse(steatosis_grade>0, "1+","0"),
         fibrosis_stage_new=ifelse(fibrosis_stage>0, "1+","0"))

md<-md%>%left_join(.,grade, by="host_subject_id")%>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  mutate(steatosis_grade_new=ifelse(condition=="NA","none",steatosis_grade_new),
         fibrosis_stage_new=ifelse(condition=="NA","none",fibrosis_stage_new))

write.table(md,"16s/stool/metadata_cln_addmoreNASHcat.txt",sep = "\t",row.names = FALSE, quote=FALSE)
##################################################################
#plot RPCA of all data

ord <- read_qza("16s/stool/rpca_results/rpca_results_NASH_all_gg2/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/stool/rpca_results/rpca_results_NASH_all_gg2/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/stool/rpca_results/rpca_results_NASH_all_gg2/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/stool/metadata_cln.txt")%>%
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

# p<-rpca %>%
#   ggplot(aes(x=PC1, y=PC2, color=condition_ZT, shape=collection_disease_stage)) +
#   geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
#   theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = condition_ZT))+
#   #scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
#   scale_shape_manual(values=c(3,16)) +
#   labs(color="condition",
#        x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
#        y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH stool 16S")+ theme(plot.title = element_text(face = "bold"))
# ggsave("16s/stool/rpca_results/rpca_results_NASH_all/SFR23_1009_NASH_stool16s_RPCA.pdf", plot=p,height=4, width=4)

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=stage_condition, shape=phase)) +
  geom_point(alpha=1.0, size=2.5) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = condition_ZT))+
  scale_color_manual(values=c("#56B4E9","#0072B2","#E69F00","#D55E00","limegreen","darkgreen"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH stool 16S")+ 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")
ggsave("16s/stool/rpca_results/rpca_results_NASH_all_gg2/SFR24_0209_NASH_stool16s_RPCA_sepwk.pdf", plot=p,height=4, width=4)

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

ggsave("16s/stool/rpca_results/rpca_results_NASH_all_gg2/SFR24_0508_NASH_stool16s_RPCA_sepwk_addbxplt.pdf", plot=final_plot,height=4.7, width=5.7)

#plotting Wk8 FAFT

ord <- read_qza("16s/stool/rpca_results/rpca_results_NASH_8wk_FAFT_gg2/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/stool/rpca_results/rpca_results_NASH_8wk_FAFT_gg2/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/stool/rpca_results/rpca_results_NASH_8wk_FAFT_gg2/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

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
ggsave("16s/stool/rpca_results/rpca_results_NASH_8wk_FAFT_gg2/SFR24_0209_NASH_stool16s_RPCA_sepLD.pdf", plot=p,height=4, width=4)


p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, fill=condition, shape=NASH_category)) +
  geom_point(alpha=1.0, size=4) + 
  theme_classic() +
  scale_fill_manual(values=c("#D55E00","#009E73"))+
  scale_shape_manual(values=c(21,24)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH stool 16S Wk8")+ 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave("16s/stool/rpca_results/rpca_results_NASH_8wk_FAFT_gg2/SFR24_0209_NASH_stool16s_RPCA_sepNASH.pdf", plot=p,height=4, width=4)

#plotting Wk12 FAFT

ord <- read_qza("16s/stool/rpca_results/rpca_results_NASH_12wk_FAFT_gg2/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/stool/rpca_results/rpca_results_NASH_12wk_FAFT_gg2/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/stool/rpca_results/rpca_results_NASH_12wk_FAFT_gg2/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

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
ggsave("16s/stool/rpca_results/rpca_results_NASH_12wk_FAFT_gg2/SFR24_0209_NASH_stool16s_RPCA_sepLD.pdf", plot=p,height=4, width=4)


p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, fill=condition, shape=NASH_category)) +
  geom_point(alpha=1.0, size=4) + 
  theme_classic() +
  scale_fill_manual(values=c("#D55E00","#009E73"))+
  scale_shape_manual(values=c(21,24)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH stool 16S Wk12")+ 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave("16s/stool/rpca_results/rpca_results_NASH_12wk_FAFT_gg2/SFR24_0209_NASH_stool16s_RPCA_sepNASH.pdf", plot=p,height=4, width=4)


#plotting ZT1 FA RPCA
ord <- read_qza("16s/stool/rpca_results/rpca_results_NASH_ZT1_FA/ordination.qza")
samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT1_FA/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT1_FA/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/stool/metadata_cln_addmoreNASHcat.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         fibrosis_stage=as.factor(fibrosis_stage),
         steatosis_grade=as.factor(steatosis_grade))

rpca<-rpca_dat(ord)
p<-rpca_plot_disease_timediff(rpca,"NASH_category", "collection_disease_stage", 2) + ggtitle("stool 16S ZT1 FA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT1_FA/SFR23_1009_NASH_stool16s_RPCA_ZT1FA.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease_timediff(rpca,"fibrosis_stage", "collection_disease_stage", 3) + ggtitle("stool 16S ZT1 FA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT1_FA/SFR23_1009_NASH_stool16s_RPCA_ZT1FA_fibrosisstg.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease_timediff(rpca,"steatosis_grade", "collection_disease_stage", 3) + ggtitle("stool 16S ZT1 FA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT1_FA/SFR23_1009_NASH_stool16s_RPCA_ZT1FA_steatosis_grd.pdf", plot=p,height=3, width=3)

##plot by disease stage
p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=collection_disease_stage, shape=collection_timepoint)) +
  geom_point(alpha=1.0, size=2.5, stroke=2) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = collection_disease_stage))+
  scale_color_manual(values=c("#E69F00","#D55E00"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("stool 16S ZT1 FA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT1_FA/SFR23_1009_NASH_stool16s_RPCA_ZT1FA_bywk.pdf", plot=p,height=4, width=4)


#plotting ZT13 FA RPCA

ord <- read_qza("16s/stool/rpca_results/rpca_results_NASH_ZT13_FA/ordination.qza")
samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT13_FA/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT13_FA/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/stool/metadata_cln_addmoreNASHcat.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         fibrosis_stage=as.factor(fibrosis_stage),
         steatosis_grade=as.factor(steatosis_grade))

rpca<-rpca_dat(ord)
p<-rpca_plot_disease_timediff(rpca,"NASH_category", "collection_disease_stage", 2) +ggtitle("stool 16S ZT13 FA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT13_FA/SFR23_1009_NASH_stool16s_RPCA_ZT13FA.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease_timediff(rpca,"fibrosis_stage", "collection_disease_stage", 3) + ggtitle("stool 16S ZT13 FA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT13_FA/SFR23_1009_NASH_stool16s_RPCA_ZT13FA_fibrosisstg.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease_timediff(rpca,"steatosis_grade", "collection_disease_stage", 3) + ggtitle("stool 16S ZT13 FA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT13_FA/SFR23_1009_NASH_stool16s_RPCA_ZT13FA_steatosis_grd.pdf", plot=p,height=3, width=3)

##plot by disease stage
p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=collection_disease_stage)) +
  geom_point(alpha=1.0, size=2.5) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = collection_disease_stage))+
  scale_color_manual(values=c("#E69F00","#D55E00"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("stool 16S ZT13 FA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT13_FA/SFR23_1009_NASH_stool16s_RPCA_ZT13FA_bywk.pdf", plot=p,height=4, width=4)


#plotting ZT1 FT RPCA

ord <- read_qza("16s/stool/rpca_results/rpca_results_NASH_ZT1_FT/ordination.qza")
samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT1_FT/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT1_FT/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/stool/metadata_cln_addmoreNASHcat.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         fibrosis_stage=as.factor(fibrosis_stage),
         steatosis_grade=as.factor(steatosis_grade))

rpca<-rpca_dat(ord)
p<-rpca_plot_disease_timediff(rpca,"NASH_category", "collection_disease_stage", 2) +ggtitle("stool 16S ZT1 FT")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT1_FT/SFR23_1009_NASH_stool16s_RPCA_ZT1FT.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease_timediff(rpca,"fibrosis_stage", "collection_disease_stage", 3) + ggtitle("stool 16S ZT1 FT")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT1_FT/SFR23_1009_NASH_stool16s_RPCA_ZT1FT_fibrosisstg.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease_timediff(rpca,"steatosis_grade", "collection_disease_stage", 3) + ggtitle("stool 16S ZT1 FT")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT1_FT/SFR23_1009_NASH_stool16s_RPCA_ZT1FT_steatosis_grd.pdf", plot=p,height=3, width=3)

##plot by disease stage
p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=collection_disease_stage, shape=collection_timepoint)) +
  geom_point(alpha=1.0, size=2.5, stroke=2) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = collection_disease_stage))+
  scale_color_manual(values=c("limegreen","darkgreen"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("stool 16S ZT1 FT")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT1_FT/SFR23_1009_NASH_stool16s_RPCA_ZT1FT_bywk.pdf", plot=p,height=4, width=4)


#plotting ZT13 FT RPCA
ord <- read_qza("16s/stool/rpca_results/rpca_results_NASH_ZT13_FT/ordination.qza")
samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT13_FT/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT13_FT/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/stool/metadata_cln_addmoreNASHcat.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         fibrosis_stage=as.factor(fibrosis_stage),
         steatosis_grade=as.factor(steatosis_grade))

rpca<-rpca_dat(ord)
p<-rpca_plot_disease_timediff(rpca,"NASH_category", "collection_disease_stage", 2) +ggtitle("stool 16S ZT13 FT")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT13_FT/SFR23_1009_NASH_stool16s_RPCA_ZT13FT.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease_timediff(rpca,"fibrosis_stage", "collection_disease_stage", 3) + ggtitle("stool 16S ZT13 FT")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT13_FT/SFR23_1009_NASH_stool16s_RPCA_ZT13FT_fibrosisstg.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease_timediff(rpca,"steatosis_grade", "collection_disease_stage", 3) + ggtitle("stool 16S ZT13 FT")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT13_FT/SFR23_1009_NASH_stool16s_RPCA_ZT13FT_steatosis_grd.pdf", plot=p,height=3, width=3)


##plot by disease stage
p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=collection_disease_stage)) +
  geom_point(alpha=1.0, size=2.5) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = collection_disease_stage))+
  scale_color_manual(values=c("limegreen","darkgreen"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("stool 16S ZT13 FT")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT13_FT/SFR23_1009_NASH_stool16s_RPCA_ZT13FT_bywk.pdf", plot=p,height=4, width=4)


#plotting ZT1 NA RPCA

ord <- read_qza("16s/stool/rpca_results/rpca_results_NASH_ZT1_NA/ordination.qza")
samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT1_NA/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT1_NA/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/stool/metadata_cln.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-rpca_dat(ord)
p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=NASH_category, shape=collection_disease_stage)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = collection_disease_stage))+
  scale_color_manual(values=c("#666666"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("stool 16S ZT1 NA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT1_NA/SFR23_1009_NASH_stool16s_RPCA_ZT1NA.pdf", plot=p,height=3, width=3)
#ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT1_NA/SFR23_1009_NASH_stool16s_RPCA_ZT1NA_welllip.pdf", plot=p,height=3, width=3)

#plotting ZT13 NA RPCA

ord <- read_qza("16s/stool/rpca_results/rpca_results_NASH_ZT13_NA/ordination.qza")
samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT13_NA/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/stool/rpca_results/rpca_results_NASH_ZT13_NA/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/stool/metadata_cln.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-rpca_dat(ord)

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=NASH_category, shape=collection_disease_stage)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = collection_disease_stage))+
  scale_color_manual(values=c("#666666"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("stool 16S ZT13 NA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT13_NA/SFR23_1009_NASH_stool16s_RPCA_ZT13NA.pdf", plot=p,height=3, width=3)
#ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT13_NA/SFR23_1009_NASH_stool16s_RPCA_ZT13NA_welllip.pdf", plot=p,height=3, width=3)

#plotting Wk12 FA RPCA

ord <- read_qza("16s/stool/rpca_results/rpca_results_NASH_12wk_FA/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/stool/rpca_results/rpca_results_NASH_12wk_FA/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/stool/rpca_results/rpca_results_NASH_12wk_FA/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("16s/stool/metadata_cln_addmoreNASHcat.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         fibrosis_stage=as.factor(fibrosis_stage),
         steatosis_grade=as.factor(steatosis_grade))

rpca<-rpca_dat(ord)
p<-rpca_plot_disease(rpca,"NASH_category",2)+ggtitle("stool 16S Wk12 FA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_12wk_FA/SFR23_1017_NASH_stool16s_RPCA_Wk12FA.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease(rpca,"fibrosis_stage",3)+ggtitle("stool 16S Wk12 FA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_12wk_FA/SFR23_1017_NASH_stool16s_RPCA_Wk12FA_fibrosisstg.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease(rpca,"steatosis_grade",3)+ggtitle("stool 16S Wk12 FA")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_12wk_FA/SFR23_1017_NASH_stool16s_RPCA_Wk12FA_steatosisgrd.pdf", plot=p,height=3, width=3)

#plotting Wk12 FT RPCA

ord <- read_qza("16s/stool/rpca_results/rpca_results_NASH_12wk_FT/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"16s/stool/rpca_results/rpca_results_NASH_12wk_FT/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"16s/stool/rpca_results/rpca_results_NASH_12wk_FT/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

rpca<-rpca_dat(ord)
p<-rpca_plot_disease(rpca,"NASH_category",2)+ggtitle("stool 16S Wk12 FT")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_12wk_FT/SFR23_1017_NASH_stool16s_RPCA_Wk12FT.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease(rpca,"fibrosis_stage",3)+ggtitle("stool 16S Wk12 FT")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_12wk_FT/SFR23_1017_NASH_stool16s_RPCA_Wk12FT_fibrosisstg.pdf", plot=p,height=3, width=3)

p<-rpca_plot_disease(rpca,"steatosis_grade",3)+ggtitle("stool 16S Wk12 FT")+ theme(plot.title = element_text(face = "bold"))
ggsave("16s/stool/rpca_results/rpca_results_NASH_12wk_FT/SFR23_1017_NASH_stool16s_RPCA_Wk12FT_steatosisgrd.pdf", plot=p,height=3, width=3)
##################################################################
md<-fread("16s/stool/metadata_cln.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(M1=sample_name)%>%
  dplyr::select(M1,condition,NASH_category,collection_timepoint,collection_disease_stage)

distance_max<-fread("16s/stool/rpca_results/rpca_results_NASH_all/distance_matrix/distance-matrix.tsv")%>%
  dplyr::rename(M1=V1)%>%gather(M2,distance,-M1)%>%
  separate(M1, c("M1_host","M1_time_point"),remove=FALSE)%>%
  separate(M2, c("M2_host","M2_time_point"),remove=FALSE)%>%
  filter(M1_host==M2_host)%>%
  filter(M1!=M2)%>%
  left_join(.,md,by="M1")%>%
  filter(collection_disease_stage=="TRF" & (M2_time_point=="N1"|M2_time_point=="N2"))%>%
  mutate(NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH")),
         condition=factor(condition,levels=c("NA","FA","FT")),
         collection_timepoint=factor(collection_timepoint, levels=c("ZT1","ZT13")))%>%
  filter(condition!="NA")

#boxplot of w/i wk8 to wk12 dist (comb ZT) to see ZT1 vs ZT13 diff (none)
p<-ggplot(distance_max, aes(x=condition, y=distance,fill=collection_timepoint)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_minimal()+ scale_color_manual(values=c("#0000a7","#c1272d"))+
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(x="condition",y="within mice robust aitchison distance", title="NASH all (wk8 to wk12 distance)")+
  theme(legend.position = "top")
ggsave("16s/stool/rpca_results/rpca_results_NASH_all/distance_matrix/NASH_all_wk8wk12dist_ZT1ZT13_wNA.pdf", plot=p,height=3.5, width=4)

perform_wilcox_test <- function(data) {pairwise.t.test(data$distance, data$collection_timepoint, p.adjust.method = "fdr")}
perform_wilcox_test2 <- function(data) {pairwise.t.test(data$distance, data$condition, p.adjust.method = "fdr")}

cond_p <- distance_max %>%
  split(.$condition) %>%
  map(perform_wilcox_test) %>%
  map_dbl(pluck, "p.value")
# NA        FA        FT 
# 0.1261699 0.6627543 0.3218857 

cond_w <- distance_max %>%
  split(.$collection_timepoint) %>%
  map(perform_wilcox_test2) 

cond_w$ZT1$p.value
# NA        FA
# FA 0.6797745        NA
# FT 0.6797745 0.6797745

cond_w$ZT13$p.value
# NA       FA
# FA 0.1478817       NA
# FT 0.4281742 0.336129

#boxplot of w/i wk8 to wk 12 dist (comb ZT) to see for non-NASH NASH diff 
#(there's some between nonNASH vsNASH for FT and FA vs FT for NASH)
p<-ggplot(distance_max, aes(x=condition, y=distance,fill=NASH_category)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+ scale_color_manual(values=c("#0000a7","#c1272d"))+
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(x="condition",y="within mice robust aitchison distance", title="NASH all (wk8 to wk12 distance)")+
  theme(legend.position = "top")
ggsave("16s/stool/rpca_results/rpca_results_NASH_all/distance_matrix/NASH_all_wk8w12dist_nonNASHNASHcat.pdf", plot=p,height=3.5, width=3.5)

perform_wilcox_test3 <- function(data) {pairwise.t.test(data$distance, data$NASH_category, p.adjust.method = "fdr")}
perform_wilcox_test4 <- function(data) {pairwise.t.test(data$distance, data$condition, p.adjust.method = "fdr")}

cond_p <- distance_max %>%
  split(.$condition) %>%
  map(perform_wilcox_test3) 

cond_p$FA$p.value
# Non_NASH
# NASH 0.1881065
cond_p$FT$p.value
# Non_NASH
# NASH 0.03514265

cond_w <- distance_max %>%
  split(.$NASH_category) %>%
  map(perform_wilcox_test4) 
cond_w$Non_NASH$p.value
# FA
# FT 0.138919
cond_w$NASH$p.value
# FA
# FT 0.05034304

#boxplot of w/i wk8 to wk12 dist (ZT1) to see for non-NASH NASH diff (there's none)

ZT1dist<-fread("16s/stool/rpca_results/rpca_results_NASH_ZT1/distance_matrix/distance-matrix.tsv")%>%
  dplyr::rename(M1=V1)%>%gather(M2,distance,-M1)%>%
  separate(M1, c("M1_host","M1_time_point"),remove=FALSE)%>%
  separate(M2, c("M2_host","M2_time_point"),remove=FALSE)%>%
  filter(M1_host==M2_host)%>%
  filter(M1!=M2)%>%
  left_join(.,md,by="M1")%>%
  filter(collection_disease_stage=="TRF")%>%
  mutate(NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH")),
         condition=factor(condition,levels=c("NA","FA","FT")))%>%
  filter(condition!="NA")

#ZT1dist<-distance_max%>%filter(collection_timepoint=="ZT1" & M2_time_point=="N1")
p<-ggplot(ZT1dist, aes(x=condition, y=distance,fill=NASH_category)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_minimal()+ scale_color_manual(values=c("#0000a7","#c1272d"))+
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(x="condition",y="within mice robust aitchison distance", title="NASH ZT1 (wk8 to wk12 distance)")+
  theme(legend.position = "top")
#ggsave("16s/stool/rpca_results/rpca_results_NASH_all/distance_matrix/NASH_ZT1_wk8w12dist_nonNASHNASHcat.pdf", plot=p,height=3.5, width=3.5)
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT1/distance_matrix/NASH_ZT1_wk8w12dist_nonNASHNASHcat.pdf", plot=p,height=3.5, width=3.5)

perform_wilcox_test3 <- function(data) {pairwise.t.test(data$distance, data$NASH_category, p.adjust.method = "fdr")}
perform_wilcox_test4 <- function(data) {pairwise.t.test(data$distance, data$condition, p.adjust.method = "fdr")}

cond_p <- ZT1dist %>%
  split(.$condition) %>%
  map(perform_wilcox_test3) 
cond_p$FA$p.value
# Non_NASH
# NASH 0.9038319
cond_p$FT$p.value
# Non_NASH
# NASH 0.1392635

cond_w <- ZT1dist %>%
  split(.$NASH_category) %>%
  map(perform_wilcox_test4) 

cond_w$Non_NASH$p.value
# FA
# FT 0.4675668
cond_w$NASH$p.value
# FA
# FT 0.2356717

#boxplot of w/i wk8 to wk12 dist (ZT12) to see for non-NASH NASH diff (there's none)

ZT13dist<-fread("16s/stool/rpca_results/rpca_results_NASH_ZT13/distance_matrix/distance-matrix.tsv")%>%
  dplyr::rename(M1=V1)%>%gather(M2,distance,-M1)%>%
  separate(M1, c("M1_host","M1_time_point"),remove=FALSE)%>%
  separate(M2, c("M2_host","M2_time_point"),remove=FALSE)%>%
  filter(M1_host==M2_host)%>%
  filter(M1!=M2)%>%
  left_join(.,md,by="M1")%>%
  filter(collection_disease_stage=="TRF")%>%
  mutate(NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH")),
         condition=factor(condition,levels=c("NA","FA","FT")))%>%
  filter(condition!="NA")

#ZT13dist<-distance_max%>%filter(collection_timepoint=="ZT13" & M2_time_point=="N2")
p<-ggplot(ZT13dist, aes(x=condition, y=distance,fill=NASH_category)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_minimal()+ scale_color_manual(values=c("#0000a7","#c1272d"))+
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(x="condition",y="within mice robust aitchison distance", title="NASH ZT13 (wk8 to wk12 distance)")+
  theme(legend.position = "top")
#ggsave("16s/stool/rpca_results/rpca_results_NASH_all/distance_matrix/NASH_ZT13_wk8w12dist_nonNASHNASHcat.pdf", plot=p,height=3.5, width=3.5)
ggsave("16s/stool/rpca_results/rpca_results_NASH_ZT13/distance_matrix/NASH_ZT13_wk8w12dist_nonNASHNASHcat.pdf", plot=p,height=3.5, width=3.5)

perform_wilcox_test3 <- function(data) {pairwise.t.test(data$distance, data$NASH_category, p.adjust.method = "fdr")}
perform_wilcox_test4 <- function(data) {pairwise.t.test(data$distance, data$condition, p.adjust.method = "fdr")}

cond_p <- ZT13dist %>%
  split(.$condition) %>%
  map(perform_wilcox_test3) 
cond_p$FA$p.value
# Non_NASH
# NASH 0.386393
cond_p$FT$p.value
# Non_NASH
# NASH 0.1581974

cond_w <- ZT13dist %>%
  split(.$NASH_category) %>%
  map(perform_wilcox_test4)
cond_w$Non_NASH$p.value
# FA
# FT 0.1213329
cond_w$NASH$p.value
# FA
# FT 0.781539
###############

md<-fread("16s/stool/metadata_cln.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(M1=sample_name)%>%
  dplyr::select(M1,condition,NASH_category,collection_timepoint,collection_disease_stage)

#by ZT time
distance_max<-fread("16s/stool/rpca_results/rpca_results_NASH_all/distance_matrix/distance-matrix.tsv")%>%
  dplyr::rename(M1=V1)%>%gather(M2,distance,-M1)%>%
  separate(M1, c("M1_host","M1_time_point"),remove=FALSE)%>%
  separate(M2, c("M2_host","M2_time_point"),remove=FALSE)%>%
  filter(M1_host==M2_host)%>%
  filter(M1!=M2)%>%
  left_join(.,md,by="M1")%>%
  filter(collection_timepoint=="ZT1" & (M2_time_point=="N2"|M2_time_point=="T2"))%>%
  mutate(NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH")),
         collection_disease_stage=factor(collection_disease_stage,levels=c("TRF","NASH")),
         condition=factor(condition,levels=c("NA","FA","FT")))%>%
  filter(condition!="NA")

#boxplot of w/i ZT1 to ZT13 dist (comb wks) to see 8 vs 12 wk diff
p<-ggplot(distance_max, aes(x=condition, y=distance,fill=collection_disease_stage)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_minimal()+ scale_color_manual(values=c("#0000a7","#c1272d"))+
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(x="condition",y="within mice robust aitchison distance", title="NASH all (ZT1 to ZT13 distance)")+
  theme(legend.position = "top")
ggsave("16s/stool/rpca_results/rpca_results_NASH_all/distance_matrix/NASH_all_ZT1ZT13dist_wk8wk12_wNA.pdf", plot=p,height=3.5, width=4)

perform_wilcox_test <- function(data) {pairwise.wilcox.test(data$distance, data$collection_disease_stage, p.adjust.method = "fdr")}
perform_wilcox_test2 <- function(data) {pairwise.wilcox.test(data$distance, data$condition, p.adjust.method = "fdr")}

cond_p <- distance_max %>%
  split(.$condition) %>%
  map(perform_wilcox_test) %>%
  map_dbl(pluck, "p.value")

# NA        FA        FT 
# 0.2189208 0.4428332 0.2415238

cond_w <- distance_max %>%
  split(.$collection_disease_stage) %>%
  map(perform_wilcox_test2) 

cond_w[["TRF"]]
# NA   FA  
# FA 0.93 -   
#   FT 0.93 0.93

cond_w[["NASH"]]
# NA FA
# FA 1  - 
#   FT 1  1 

#boxplot of w/i ZT1 to ZT13 dist (comb wks) to see for non-NASH NASH diff (there's none)
p<-ggplot(distance_max, aes(x=condition, y=distance,fill=NASH_category)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_minimal()+ scale_color_manual(values=c("#0000a7","#c1272d"))+
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(x="condition",y="within mice robust aitchison distance", title="NASH all (ZT1 to ZT13 distance)")+
  theme(legend.position = "top")
ggsave("16s/stool/rpca_results/rpca_results_NASH_all/distance_matrix/NASH_all_ZT1ZT13dist_nonNASHNASHcat.pdf", plot=p,height=3.5, width=3.5)

perform_wilcox_test3 <- function(data) {pairwise.wilcox.test(data$distance, data$NASH_category, p.adjust.method = "fdr")}
perform_wilcox_test4 <- function(data) {pairwise.wilcox.test(data$distance, data$condition, p.adjust.method = "fdr")}

cond_p <- distance_max %>%
  split(.$condition) %>%
  map(perform_wilcox_test3)
cond_p$FA$p.value
# Non_NASH
# NASH 0.2381413
cond_p$FT$p.value
# Non_NASH
# NASH 0.2635644

cond_w <- distance_max %>%
  split(.$NASH_category) %>%
  map(perform_wilcox_test4) 
cond_w$Non_NASH$p.value
# FA
# FT 0.2344988
cond_w$NASH$p.value
# FA
# FT 0.2240065

#boxplot of w/i ZT1 to ZT13 dist (wk8) to see for non-NASH NASH diff (there's none)
distance_max_sub<-distance_max%>%filter(collection_disease_stage=="TRF"& M2_time_point=="T2")

p<-ggplot(distance_max_sub, aes(x=condition, y=distance,fill=NASH_category)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_minimal()+ scale_color_manual(values=c("#0000a7","#c1272d"))+
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(x="condition",y="within mice robust aitchison distance", title="NASH wk8 (ZT1 to ZT13 distance)")+
  theme(legend.position = "top")
ggsave("16s/stool/rpca_results/rpca_results_NASH_all/distance_matrix/NASH_8wk_ZT1ZT13dist_nonNASHNASHcat.pdf", plot=p,height=3.5, width=3.5)

cond_p <- distance_max_sub %>%
  split(.$condition) %>%
  map(perform_wilcox_test3)
cond_p$FA$p.value
# Non_NASH
# NASH 0.1333333
cond_p$FT$p.value
# Non_NASH
# NASH 0.5333333

cond_w <- distance_max_sub %>%
  split(.$NASH_category) %>%
  map(perform_wilcox_test4) 
cond_w$Non_NASH$p.value
# FA
# FT 0.3333333
cond_w$NASH$p.value
# FA
# FT 0.1142857

#boxplot of w/i ZT1 to ZT13 dist (wk12) to see for non-NASH NASH diff (there's none)
distance_max_sub<-distance_max%>%filter(collection_disease_stage=="NASH"& M2_time_point=="N2")

p<-ggplot(distance_max_sub, aes(x=condition, y=distance,fill=NASH_category)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_minimal()+ scale_color_manual(values=c("#0000a7","#c1272d"))+
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(x="condition",y="within mice robust aitchison distance", title="NASH wk8 (ZT1 to ZT13 distance)")+
  theme(legend.position = "top")
ggsave("16s/stool/rpca_results/rpca_results_NASH_all/distance_matrix/NASH_12wk_ZT1ZT13dist_nonNASHNASHcat.pdf", plot=p,height=3.5, width=3.5)

cond_p <- distance_max_sub %>%
  split(.$condition) %>%
  map(perform_wilcox_test3)

cond_p$FA$p.value
# Non_NASH
# NASH        1
cond_p$FT$p.value
# Non_NASH
# NASH 0.2666667

cond_w <- distance_max_sub %>%
  split(.$NASH_category) %>%
  map(perform_wilcox_test4) 

cond_w$Non_NASH$p.value
# FA
# FT 0.3333333
cond_w$NASH$p.value
# FA
# FT 0.1142857

##################################################################

#plot BIRDMAn results

##non-NASH VS. NASH wk8 FT (w/i distance calc showed diff)
# annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
#   dplyr::rename(FeatureID=Feature.ID)

annot<-read_qza("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/04.sklearn.gg2.asv.taxonomy.qza")$data%>%
  dplyr::rename(FeatureID=Feature.ID)
write.table(annot,"files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/taxonomy.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

FT<-fread("16s/stool/birdman_results/NASH_TRF_FT_taxonomy_filtered.gg2.asv.counts.beta_var.tsv")
FT<-birdman_dat(FT,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FT,"16s/stool/birdman_results/birdman_justcred_wk8FT_nonNASHNASH_results_gg2.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FT) + labs(title="FT Wk8\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk8FT_nonNASHNASH_gg2.pdf",height=6, width=13)

FT$Name <- factor(FT$Name,levels = FT$Name)
ggplot(FT, aes(x =Name , y = ratio, ymin = min, ymax = max, fill=cat)) + 
  geom_bar(stat="identity", alpha=0.7)+
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  #geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(title="FT Wk8\nnon-NASH vs. NASH")+
  coord_flip()
ggsave("16s/stool/birdman_results/birdman_justcred_Wk8FT_nonNASHNASH_colored_gg2.pdf",height=6, width=14)

##non-NASH VS. NASH wk12 FT (w/i distance calc showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/taxonomy.tsv")

FT<-fread("16s/stool/birdman_results/NASH_NASH_FT_taxonomy_filtered.gg2.asv.counts.beta_var.tsv")
FT<-birdman_dat(FT,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FT,"16s/stool/birdman_results/birdman_justcred_wk12FT_nonNASHNASH_results_gg2.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FT) + labs(title="FT Wk12\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk12FT_nonNASHNASH_gg2.pdf",height=4, width=13)

FT$Name <- factor(FT$Name,levels = FT$Name)
ggplot(FT, aes(x =Name , y = ratio, ymin = min, ymax = max, fill=cat)) + 
  geom_bar(stat="identity", alpha=0.7)+
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  #geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(title="FT Wk12\nnon-NASH vs. NASH")+
  coord_flip()
ggsave("16s/stool/birdman_results/birdman_justcred_Wk12FT_nonNASHNASH_colored_gg2.pdf",height=4, width=14)

#create an alluvial plot of the FT data for non-nash vs nash going from wk 8 to wk 12

FTw12<-fread("16s/stool/birdman_results/birdman_justcred_wk12FT_nonNASHNASH_results_gg2.tsv") %>%
  arrange(desc(ratio))%>%rownames_to_column("rank")%>%mutate(week="Wk12",
                                                             Name_rank=paste(rank,Name, sep="."))

FTw8<-fread("16s/stool/birdman_results/birdman_justcred_wk8FT_nonNASHNASH_results_gg2.tsv")%>%
  arrange(desc(ratio))%>%rownames_to_column("rank")%>%mutate(week="Wk8",
                                                             Name_rank=paste(rank,Name, sep="."))

FT<-rbind(FTw8,FTw12)%>%
  mutate(week=factor(week, levels=c("Wk8","Wk12")),
         cat=factor(cat,levels=c("Non_NASH","NASH")),
         rank=as.numeric(rank))

FT$Name=fct_reorder(FT$Name, FT$rank, max)
FT$Name_rank=fct_reorder(FT$Name_rank, FT$rank, max)

ggplot(FT, aes(x = week,fill = cat, stratum = Name_rank, alluvium = Name)) +
  stat_stratum(alpha = .25) +
  geom_flow()+
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  stat_alluvium(geom = "text", size=3,aes(label = Name_rank))+theme_minimal()

ggsave("16s/stool/birdman_results/birdman_justcred_wk8wk12FT_alluvial_gg2.pdf",height=4, width=8)


##non-NASH VS. NASH wk8 FA (w/i distance calc showed no diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/taxonomy.tsv")

FA<-fread("16s/stool/birdman_results/NASH_TRF_FA_taxonomy_filtered.gg2.asv.counts.beta_var.tsv")
FA<-birdman_dat(FA,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FA,"16s/stool/birdman_results/birdman_justcred_wk8FA_nonNASHNASH_results_gg2.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FA) + labs(title="FA Wk8\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk8FA_nonNASHNASH_gg2.pdf",height=6, width=14)

FA$Name <- factor(FA$Name,levels = FA$Name)
ggplot(FA, aes(x =Name , y = ratio, ymin = min, ymax = max, fill=cat)) + 
  geom_bar(stat="identity", alpha=0.7)+
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  #geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(title="FA Wk8\nnon-NASH vs. NASH")+
  coord_flip()
ggsave("16s/stool/birdman_results/birdman_justcred_Wk8FA_nonNASHNASH_colored_gg2.pdf",height=6, width=15)

##non-NASH VS. NASH wk12 FA (w/i distance calc showed no diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/taxonomy.tsv")

FA<-fread("16s/stool/birdman_results/NASH_NASH_FA_taxonomy_filtered.gg2.asv.counts.beta_var.tsv")
FA<-birdman_dat(FA,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FA,"16s/stool/birdman_results/birdman_justcred_wk12FA_nonNASHNASH_results_gg2.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FA) + labs(title="FA Wk12\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk12FA_nonNASHNASH_gg2.pdf",height=7, width=14)

FA$Name <- factor(FA$Name,levels = FA$Name)
ggplot(FA, aes(x =Name , y = ratio, ymin = min, ymax = max, fill=cat)) + 
  geom_bar(stat="identity", alpha=0.7)+
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  #geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(title="FA Wk12\nnon-NASH vs. NASH")+
  coord_flip()
ggsave("16s/stool/birdman_results/birdman_justcred_Wk12FA_nonNASHNASH_colored_gg2.pdf",height=7, width=15)

#create an alluvial plot of the FA data for non-nash vs nash going from wk 8 to wk 12

FAw12<-fread("16s/stool/birdman_results/birdman_justcred_wk12FA_nonNASHNASH_results_gg2.tsv") %>%
  arrange(desc(ratio))%>%rownames_to_column("rank")%>%mutate(week="Wk12",
                                                             Name_rank=paste(rank,Name, sep="."))

FAw8<-fread("16s/stool/birdman_results/birdman_justcred_wk8FA_nonNASHNASH_results_gg2.tsv")%>%
  arrange(desc(ratio))%>%rownames_to_column("rank")%>%mutate(week="Wk8",
                                                             Name_rank=paste(rank,Name, sep="."))

FA<-rbind(FAw8,FAw12)%>%
  mutate(week=factor(week, levels=c("Wk8","Wk12")),
         cat=factor(cat,levels=c("Non_NASH","NASH")),
         rank=as.numeric(rank))

FA$Name=fct_reorder(FA$Name, FA$rank, max)
FA$Name_rank=fct_reorder(FA$Name_rank, FA$rank, max)

ggplot(FA, aes(x = week,fill = cat, stratum = Name_rank, alluvium = Name)) +
  stat_stratum(alpha = .25) +
  geom_flow()+
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  stat_alluvium(geom = "text", size=3,aes(label = Name_rank))+theme_minimal()

ggsave("16s/stool/birdman_results/birdman_justcred_wk8wk12FA_alluvial_gg2.pdf",height=5, width=8)

#####################

##non-NASH VS. NASH wk8 FA ZT1 (w/i distance calc showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

FAZ1<-fread("16s/stool/birdman_results/NASH_TRF_FA_ZT1_taxonomy_filtered.asv.counts.beta_var.tsv")
FAZ1<-birdman_dat(FAZ1,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FAZ1,"16s/stool/birdman_results/birdman_justcred_wk8FAZT1_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FAZ1) + labs(title="FA Wk8 ZT1\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk8FAZT1_nonNASHNASH.pdf",height=3, width=11)

##non-NASH VS. NASH wk8 FA ZT13 (w/i distance calc showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

FAZ13<-fread("16s/stool/birdman_results/NASH_TRF_FA_ZT13_taxonomy_filtered.asv.counts.beta_var.tsv")
FAZ13<-birdman_dat(FAZ13,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FAZ13,"16s/stool/birdman_results/birdman_justcred_wk8FAZT13_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FAZ13) + labs(title="FA Wk8 ZT13\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk8FAZT13_nonNASHNASH.pdf",height=4, width=11)

##non-NASH VS. NASH wk8 FT ZT1 (w/i distance calc showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

FTZ1<-fread("16s/stool/birdman_results/NASH_TRF_FT_ZT1_taxonomy_filtered.asv.counts.beta_var.tsv")
FTZ1<-birdman_dat(FTZ1,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FTZ1,"16s/stool/birdman_results/birdman_justcred_wk8FTZT1_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FTZ1) + labs(title="FT Wk8 ZT1\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk8FTZT1_nonNASHNASH.pdf",height=4, width=11)

##non-NASH VS. NASH wk8 FT ZT13 (w/i distance calc showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

FTZ13<-fread("16s/stool/birdman_results/NASH_TRF_FT_ZT13_taxonomy_filtered.asv.counts.beta_var.tsv")
FTZ13<-birdman_dat(FTZ13,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FTZ13,"16s/stool/birdman_results/birdman_justcred_wk8FTZT13_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FTZ13) + labs(title="FT Wk8 ZT13\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk8FTZT13_nonNASHNASH.pdf",height=2, width=11)

##non-NASH VS. NASH wk12 FA ZT1 (w/i distance calc showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

FAZ1<-fread("16s/stool/birdman_results/NASH_NASH_FA_ZT1_taxonomy_filtered.asv.counts.beta_var.tsv")
FAZ1<-birdman_dat(FAZ1,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FAZ1,"16s/stool/birdman_results/birdman_justcred_wk12FAZT1_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FAZ1) + labs(title="FA Wk12 ZT1\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk8FAZT1_nonNASHNASH.pdf",height=5, width=11)

##non-NASH VS. NASH wk12 FA ZT13 (w/i distance calc showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

FAZ13<-fread("16s/stool/birdman_results/NASH_NASH_FA_ZT13_taxonomy_filtered.asv.counts.beta_var.tsv")
FAZ13<-birdman_dat(FAZ13,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FAZ13,"16s/stool/birdman_results/birdman_justcred_wk12FAZT13_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FAZ13) + labs(title="FA Wk12 ZT13\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk12FAZT13_nonNASHNASH.pdf",height=5, width=11)

##non-NASH VS. NASH wk8 FT ZT1 (w/i distance calc showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

FTZ1<-fread("16s/stool/birdman_results/NASH_NASH_FT_ZT1_taxonomy_filtered.asv.counts.beta_var.tsv")
FTZ1<-birdman_dat(FTZ1,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FTZ1,"16s/stool/birdman_results/birdman_justcred_wk12FTZT1_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FTZ1) + labs(title="FT Wk12 ZT1\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk12FTZT1_nonNASHNASH.pdf",height=3, width=11)

##non-NASH VS. NASH wk8 FT ZT13 (w/i distance calc showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

FTZ13<-fread("16s/stool/birdman_results/NASH_NASH_FT_ZT13_taxonomy_filtered.asv.counts.beta_var.tsv")
FTZ13<-birdman_dat(FTZ13,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FTZ13,"16s/stool/birdman_results/birdman_justcred_wk12FTZT13_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FTZ13) + labs(title="FT Wk12 ZT13\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_Wk12FTZT13_nonNASHNASH.pdf",height=4, width=11)

###################

##NASH FA VS. NASH FT (w/i distance calc showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

NASHct<-fread("16s/stool/birdman_results/NASH_NASH_FAFT_NASHcat_taxonomy_filtered.asv.counts.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(condition, Treatment('FA'))[T.FT]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(condition, Treatment('FA'))[T.FT]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('FA'))[T.FT]_hdi`))%>%
  separate(`C(condition, Treatment('FA'))[T.FT]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  left_join(.,annot, by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon, sep=" "))%>%
  filter(cred=="credible")%>%
  dplyr::select(FeatureID, ratio, min, max, Name)%>%
  arrange(ratio)

write.table(NASHct,"16s/stool/birdman_results/birdman_justcred_Wk12NASH_FAFT_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

NASHct$Name <- factor(NASHct$Name,levels = NASHct$Name)

ggplot(NASHct, aes(x =Name , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title="Wk12 NASH FA vs. FT")+
  coord_flip()

ggsave("16s/stool/birdman_results/birdman_justcred_Wk12NASH_FAFT.pdf",height=8, width=14)

###################

##non-NASH VS. NASH ZT1 FA (RPCA showed no diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

ZT1FA<-fread("16s/stool/birdman_results/NASH_ZT1_FA_taxonomy_filtered.asv.counts.beta_var.tsv")
ZT1FA<-birdman_dat(ZT1FA,annot)%>%
  mutate(ZT_time="ZT1")
write.table(ZT1FA,"16s/stool/birdman_results/birdman_justcred_ZT1FA_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT1FA) +   labs(title="FA ZT1\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_ZT1FA_nonNASHNASH.pdf",height=6, width=13.5)

##non-NASH VS. NASH ZT13 FA (RPCA showed diff)
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/results/taxonomy.csv")%>%
  dplyr::rename(FeatureID=Feature.ID)

ZT13FA<-fread("16s/stool/birdman_results/NASH_ZT13_FA_taxonomy_filtered.asv.counts.beta_var.tsv")
ZT13FA<-birdman_dat(ZT13FA,annot)%>%
  mutate(ZT_time="ZT13")
write.table(ZT13FA,"16s/stool/birdman_results/birdman_justcred_ZT13FA_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT13FA) +   labs(title="FA ZT13\nnon-NASH vs. NASH")
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

ZT1FT<-fread("16s/stool/birdman_results/NASH_ZT1_FT_taxonomy_filtered.asv.counts.beta_var.tsv")
ZT1FT<-birdman_dat(ZT1FT,annot)%>%
  mutate(ZT_time="ZT1")
write.table(ZT1FT,"16s/stool/birdman_results/birdman_justcred_ZT1FT_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT1FT) +   labs(title="FA ZT13\nnon-NASH vs. NASH")
ggsave("16s/stool/birdman_results/birdman_justcred_ZT13FA_nonNASHNASH.pdf",height=6, width=14.5)


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

ZT13FT<-fread("16s/stool/birdman_results/NASH_ZT13_FT_taxonomy_filtered.asv.counts.beta_var.tsv")
ZT13FT<-birdman_dat(ZT13FT,annot)%>%
  mutate(ZT_time="ZT13")
write.table(ZT13FT,"16s/stool/birdman_results/birdman_justcred_ZT13FT_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT13FT) +   labs(title="FT ZT13\nnon-NASH vs. NASH")
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
  arrange(ratio)%>%
  mutate(cat=ifelse(ratio<0,"Wk8","Wk12"))%>%
  mutate(cat=factor(cat,levels=c("Wk8","Wk12")))

write.table(ZT13FT,"16s/stool/birdman_results/birdman_justcred_ZT13FT_8wkv12wk_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

ZT13FT$Name <- factor(ZT13FT$Name,levels = ZT13FT$Name)

ggplot(ZT13FT, aes(x =Name , y = ratio, ymin = min, ymax = max, fill=cat)) + 
  geom_bar(stat="identity", alpha=0.7)+
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(values=c("limegreen","darkgreen"))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title="FT ZT13\n8week vs. 12 week")+
  coord_flip()

ggsave("16s/stool/birdman_results/birdman_justcred_ZT13FT_8wkv12wk.pdf",height=6, width=14.5)


#subset plot further

ZT13FT_sub<-ZT13FT[c(1:7,29:31),]
ZT13FT_sub$Name <- factor(ZT13FT_sub$Name,levels = ZT13FT_sub$Name)


ggplot(ZT13FT_sub, aes(x =Name , y = ratio, ymin = min, ymax = max, fill=cat)) +
  geom_bar(stat="identity", alpha=0.7)+
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(values=c("limegreen","darkgreen"))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title="FT ZT13\n8week vs. 12 week")+
  coord_flip()
  
  
ggsave("16s/stool/birdman_results/birdman_justcred_select10_ZT13FT_8wkv12wk.pdf",height=4, width=14.5)

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

##ZT1 vs ZT13 for FT

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

#comparing all ZT1 and ZT13 to FA and FT
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

##Wk8 FA vs FT 

wk8FA<-fread("16s/stool/birdman_results/birdman_justcred_wk8FA_nonNASHNASH_results_gg2.tsv")
wk8FT<-fread("16s/stool/birdman_results/birdman_justcred_wk8FT_nonNASHNASH_results_gg2.tsv")

list_venn <- list(FA = wk8FA$Name,
                  FT = wk8FT$Name)

p<-venn.diagram(list_venn, fill = c("#D55E00","#009E73"),height = 10,lty = 0, 
                main="NASH Wk8",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_nonNASHNASH_wk8FAFTresults_gg2.pdf")
grid.draw(p)
dev.off()


##Wk12 FA vs FT 

wk12FA<-fread("16s/stool/birdman_results/birdman_justcred_wk12FA_nonNASHNASH_results_gg2.tsv")
wk12FT<-fread("16s/stool/birdman_results/birdman_justcred_wk12FT_nonNASHNASH_results_gg2.tsv")

list_venn <- list(FA = wk12FA$Name,
                  FT = wk12FT$Name)

p<-venn.diagram(list_venn, fill = c("#D55E00","#009E73"),height = 10,lty = 0, 
                main="NASH Wk12",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_nonNASHNASH_wk12FAFTresults_gg2.pdf")
grid.draw(p)
dev.off()


##FA Wk8 to Wk12

wk8FA<-fread("16s/stool/birdman_results/birdman_justcred_wk8FA_nonNASHNASH_results_gg2.tsv")
wk12FA<-fread("16s/stool/birdman_results/birdman_justcred_wk12FA_nonNASHNASH_results_gg2.tsv")

list_venn <- list(Wk8 = wk8FA$Name,
                  Wk12 = wk12FA$Name)

p<-venn.diagram(list_venn, fill = c("#E69F00","#D55E00"),height = 10,lty = 0, 
                main="FA Non-NASH vs. NASH",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_nonNASHNASH_FAWk8Wk12results_gg2.pdf")
grid.draw(p)
dev.off()


##FT Wk8 to Wk12

wk8FT<-fread("16s/stool/birdman_results/birdman_justcred_wk8FT_nonNASHNASH_results_gg2.tsv")
wk12FT<-fread("16s/stool/birdman_results/birdman_justcred_wk12FT_nonNASHNASH_results_gg2.tsv")

list_venn <- list(Wk8 = wk8FT$Name,
                  Wk12 = wk12FT$Name)

p<-venn.diagram(list_venn, fill = c("limegreen","darkgreen"),height = 10,lty = 0, 
                main="FT non-NASH vs. NASH",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_nonNASHNASH_FTWk8Wk12results_gg2.pdf")
grid.draw(p)
dev.off()

###################
#plotting just the overlap microvbes
wk8FA<-fread("16s/stool/birdman_results/birdman_justcred_wk8FA_nonNASHNASH_results_gg2.tsv")
wk12FA<-fread("16s/stool/birdman_results/birdman_justcred_wk12FA_nonNASHNASH_results_gg2.tsv")

list_venn <- list(Wk8 = wk8FA$FeatureID,
                  Wk12 = wk12FA$FeatureID)

wk8FT<-fread("16s/stool/birdman_results/birdman_justcred_wk8FT_nonNASHNASH_results_gg2.tsv")
wk12FT<-fread("16s/stool/birdman_results/birdman_justcred_wk12FT_nonNASHNASH_results_gg2.tsv")

list_venn <- list(Wk8 = wk8FT$FeatureID,
                  Wk12 = wk12FT$FeatureID)

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

FAsel<-all$`Wk8:Wk12`
FTsel<-all$`Wk8:Wk12`

md<-fread("16s/stool/metadata_cln_addmoreNASHcat.txt")
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/taxonomy.tsv")

m16s<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/NASH_taxonomy_filtered.gg2.asv.counts-rclr.txt")%>%
  dplyr::rename(FeatureID=`#OTU ID`)%>%
  gather(sample_name,counts,-FeatureID)%>%
  mutate(counts=ifelse(is.na(counts),0,counts))%>%
  left_join(.,md,by="sample_name")%>%
  left_join(.,annot,by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon,sep=" "))

write.table(m16s,"16s/stool/quant_table_NASH_gg2_rclr_wmd.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtb_sel1<-m16s%>%filter(FeatureID %in% FAsel)%>%
  filter(collection_disease_stage=="TRF" & condition=="FA")

mtb_sel2<-m16s%>%filter(FeatureID %in% FTsel)%>%
  filter(collection_disease_stage=="TRF" & condition=="FT")

mtb_sel<-rbind(mtb_sel1,mtb_sel2)%>%
  mutate(Name=factor(Name,levels=c("65bc2b52d9d09bde0d63aadfcffa450f d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Lachnospirales; f__Lachnospiraceae",
                                   "726fe9527294f50f4d2654b3a25730a9 d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Lachnospirales; f__Lachnospiraceae; g__Ventrisoma; s__Ventrisoma faecale",
                                   "7f782768e8d3067af0a923af092239d3 d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Lachnospirales; f__Lachnospiraceae",
                                   "eb4104bfe7ff965e7288b9e381317aa8 d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Lachnospirales; f__Lachnospiraceae",
                                   "e613626b2ce8bfbf28995d9b11edaeed d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Lachnospirales; f__Lachnospiraceae; g__; s__",
                                   "d90a306bb67d75602d9e789894d2c81e d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Oscillospirales; f__Oscillospiraceae_88309; g__Lawsonibacter",
                                   "d1c1529acbbd977e6219701a664424f0 d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Lachnospirales; f__Lachnospiraceae; g__Sporofaciens; s__",
                                   "923a477c342aa445c4eecca0f8415ed5 d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Lachnospirales; f__Lachnospiraceae; g__Merdisoma; s__Merdisoma sp011959465",
                                   "a6c1d5c393eaa8b78baddc67c010b5fe d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Oscillospirales; f__Oscillospiraceae_88309; g__; s__",
                                   "be65daef90f69cfc35f5bcd8e387d377 d__Bacteria; p__Firmicutes_A; c__Clostridia_258483; o__Lachnospirales; f__Lachnospiraceae; g__Acetatifactor; s__Acetatifactor sp011959105")))

p <- ggplot(mtb_sel, aes(x=Name, y=log10(counts+1), fill=NASH_category)) + 
  geom_violin(trim=FALSE,alpha=0.7)  + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), size=0.5)+
  theme_classic() + 
  facet_grid(~condition,scales="free") +
  force_panelsizes(cols = c(0.8, 0.4)) +
  theme(legend.position = "top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("#0000a7","#c1272d","grey90")) + 
  labs(title="Wk8",x="NASH Category",y="rclr")

ggsave("16s/stool/birdman_results/birdman_shared11hits_wk8_violinplt.pdf", plot=p,height=12, width=9)
ggsave("16s/stool/birdman_results/birdman_shared11hits_wk8_violinplt_log10.pdf", plot=p,height=12, width=9)
ggsave("16s/stool/birdman_results/birdman_shared11hits_wk8_violinplt_orgw12.pdf", plot=p,height=11.5, width=10)


mtb_sel1<-m16s%>%filter(FeatureID %in% FAsel)%>%
  filter(collection_disease_stage=="NASH" & condition=="FA")

mtb_sel2<-m16s%>%filter(FeatureID %in% FTsel)%>%
  filter(collection_disease_stage=="NASH" & condition=="FT")

mtb_sel<-rbind(mtb_sel1,mtb_sel2)%>%
  mutate(Name=factor(Name,levels=c("65bc2b52d9d09bde0d63aadfcffa450f k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae",
                                   "7f782768e8d3067af0a923af092239d3 k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__; s__",
                                   "726fe9527294f50f4d2654b3a25730a9 k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae",
                                   "eb4104bfe7ff965e7288b9e381317aa8 k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Dorea; s__",
                                   "d1c1529acbbd977e6219701a664424f0 k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__[Ruminococcus]; s__gnavus",
                                   "e613626b2ce8bfbf28995d9b11edaeed k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae",
                                   "d90a306bb67d75602d9e789894d2c81e k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Oscillospira; s__",
                                   "a6c1d5c393eaa8b78baddc67c010b5fe k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Oscillospira; s__",
                                   "be65daef90f69cfc35f5bcd8e387d377 k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae",
                                   "923a477c342aa445c4eecca0f8415ed5 k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Coprococcus; s__")))

p <- ggplot(mtb_sel, aes(x=Name, y=counts, fill=NASH_category)) + 
  #geom_boxplot(alpha=0.7)+
  geom_violin(trim=FALSE,alpha=0.7)  + 
  #geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), size=0.5)+
  theme_classic() + 
  facet_grid(~condition,scales="free") +
  force_panelsizes(cols = c(0.8, 0.4)) +
  theme(legend.position = "top",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("#0000a7","#c1272d","grey90")) + 
  labs(title="Wk12",x="NASH Category",y="rclr")
ggsave("16s/stool/birdman_results/birdman_shared11hits_wk12_violinplt.pdf", plot=p,height=11.5, width=10)
ggsave("16s/stool/birdman_results/birdman_shared11hits_wk12_bxplt.pdf", plot=p,height=12, width=10)
###################
#showing paired diff for select microbes

md<-fread("16s/stool/metadata_cln_addmoreNASHcat.txt")
annot<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/taxonomy.tsv")
m16s<-fread("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/16S/stool/data/rrichter_preprocessed_20211020_ID_13785_gg2/NASH_taxonomy_filtered.gg2.asv.counts-rclr.txt")%>%
  dplyr::rename(FeatureID=`#OTU ID`)%>%
  gather(sample_name,counts,-FeatureID)%>%
  mutate(counts=ifelse(is.na(counts),0,counts))%>%
  left_join(.,md,by="sample_name")%>%
  left_join(.,annot,by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon,sep=" "),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH")),
         collection_disease_stage=factor(collection_disease_stage,levels=c("TRF","NASH")))

micro1<-paired_plots_ind(m16s,"65bc2b52d9d09bde0d63aadfcffa450f")
ggsave("16s/stool/birdman_results/bdmFA7_65bc2_nonNASHNASH_gg2_unprd.pdf",height=4 , width=4)

micro2<-paired_plots_ind(m16s,"7f782768e8d3067af0a923af092239d3")
ggsave("16s/stool/birdman_results/bdmFA7_7f782_nonNASHNASH_gg2_unprd.pdf",height=4 , width=4)

micro3<-paired_plots_ind(m16s,"726fe9527294f50f4d2654b3a25730a9")
ggsave("16s/stool/birdman_results/bdmFA7_726fe_nonNASHNASH_gg2_unprd.pdf",height=4 , width=4)

micro4<-paired_plots_ind(m16s,"eb4104bfe7ff965e7288b9e381317aa8")
ggsave("16s/stool/birdman_results/bdmFA7_eb410_nonNASHNASH_gg2_unprd.pdf",height=4 , width=4)

micro5<-paired_plots_ind(m16s,"d1c1529acbbd977e6219701a664424f0")
ggsave("16s/stool/birdman_results/bdmFA7_d1c15_nonNASHNASH_gg2_unprd.pdf",height=4 , width=4)

micro6<-paired_plots_ind(m16s,"e613626b2ce8bfbf28995d9b11edaeed")
ggsave("16s/stool/birdman_results/bdmFA7_e6136_nonNASHNASH_gg2_unprd.pdf",height=4 , width=4)

micro7<-paired_plots_ind(m16s,"d90a306bb67d75602d9e789894d2c81e")
ggsave("16s/stool/birdman_results/bdmFA7_d90a3_nonNASHNASH_gg2_unprd.pdf",height=4 , width=4)



micro1<-paired_plots_ind(m16s,"d1c1529acbbd977e6219701a664424f0")
ggsave("16s/stool/birdman_results/bdmFT4_d1c15_nonNASHNASH_gg2_unprd.pdf",height=4 , width=4)

micro2<-paired_plots_ind(m16s,"a6c1d5c393eaa8b78baddc67c010b5fe")
ggsave("16s/stool/birdman_results/bdmFT4_a6c1d_nonNASHNASH_gg2_unprd.pdf",height=4 , width=4)

micro3<-paired_plots_ind(m16s,"be65daef90f69cfc35f5bcd8e387d377")
ggsave("16s/stool/birdman_results/bdmFT4_be65d_nonNASHNASH_gg2_unprd.pdf",height=4 , width=4)

micro4<-paired_plots_ind(m16s,"923a477c342aa445c4eecca0f8415ed5")
ggsave("16s/stool/birdman_results/bdmFT4_923a4_nonNASHNASH_gg2_unprd.pdf",height=4 , width=4)

#are within w8 w12 diff significant



###################
#comparing all Wk8 and wk12 to FA and FT
list_venn <- list(Wk8_FA = wk8FA$Name,
                  Wk12_FA = wk12FA$Name,
                  Wk8_FT = wk8FT$Name,
                  Wk12_FT = wk12FT$Name)

p<-venn.diagram(list_venn,height = 10,
                main="NASH",
                width = 10,alpha = c(0.5, 0.5, 0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_nonNASHNASH_wk8wk12_results.pdf")
grid.draw(p)
dev.off()

##ZT1 vs ZT13 for FA Wk8

ZT1FA<-fread("16s/stool/birdman_results/birdman_justcred_wk8FAZT1_nonNASHNASH_results.tsv")
ZT13FA<-fread("16s/stool/birdman_results/birdman_justcred_wk8FAZT13_nonNASHNASH_results.tsv")

p<-venn_plt(ZT1FA,ZT13FA,"NASH Wk8 FA")
grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_FAwk8_nonNASHNASH_results.pdf")
grid.draw(p)
dev.off()


##ZT1 vs ZT13 for FT Wk8

ZT1FT<-fread("16s/stool/birdman_results/birdman_justcred_wk8FTZT1_nonNASHNASH_results.tsv")
ZT13FT<-fread("16s/stool/birdman_results/birdman_justcred_wk8FTZT13_nonNASHNASH_results.tsv")

p<-venn_plt(ZT1FT,ZT13FT,"NASH Wk 8 FT")
grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_FTwk8_nonNASHNASH_results.pdf")
grid.draw(p)
dev.off()

#comparing all ZT1 and ZT13 to FA and FT
list_venn <- list(ZT1_FA = ZT1FA$Name,
                  ZT13_FA = ZT13FA$Name,
                  ZT1_FT = ZT1FT$Name,
                  ZT13_FT = ZT13FT$Name)

p<-venn.diagram(list_venn,height = 10,
                main="NASH Wk8",
                width = 10,alpha = c(0.5, 0.5, 0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_nonNASHNASH_wk8_ZT1ZT13_results.pdf")
grid.draw(p)
dev.off()

##ZT1 vs ZT13 for FA Wk12

ZT1FA<-fread("16s/stool/birdman_results/birdman_justcred_wk12FAZT1_nonNASHNASH_results.tsv")
ZT13FA<-fread("16s/stool/birdman_results/birdman_justcred_wk12FAZT13_nonNASHNASH_results.tsv")

p<-venn_plt(ZT1FA,ZT13FA,"NASH Wk 12 FA")
grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_FAwk12_nonNASHNASH_results.pdf")
grid.draw(p)
dev.off()

##ZT1 vs ZT13 for FA Wk12

ZT1FT<-fread("16s/stool/birdman_results/birdman_justcred_wk12FTZT1_nonNASHNASH_results.tsv")
ZT13FT<-fread("16s/stool/birdman_results/birdman_justcred_wk12FTZT13_nonNASHNASH_results.tsv")

p<-venn_plt(ZT1FT,ZT13FT,"NASH Wk 12 FT")
grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_FTwk12_nonNASHNASH_results.pdf")
grid.draw(p)
dev.off()

#comparing all ZT1 and ZT13 to FA and FT
list_venn <- list(ZT1_FA = ZT1FA$Name,
                  ZT13_FA = ZT13FA$Name,
                  ZT1_FT = ZT1FT$Name,
                  ZT13_FT = ZT13FT$Name)

p<-venn.diagram(list_venn,height = 10,
                main="NASH Wk12",
                width = 10,alpha = c(0.5, 0.5, 0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="16s/stool/birdman_results/birdman_justcredoverlap_nonNASHNASH_wk12_ZT1ZT13_results.pdf")
grid.draw(p)
dev.off()







