setwd("~/Notebooks/sfloresr/MASH-TRF/MASHomics/")

library(tidyverse)
library(data.table)
library(ggpubr)
library(VennDiagram)
library("qiime2R")

##################################################################
#inputs
mash_metadata<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/metadata_cln_addmoreNASHcat.txt"
taxonomy<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/mashstool16s_preprocessed_20211020_ID_13785_gg2/04.sklearn.gg2.asv.taxonomy.qza"
feattab<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/mashstool16s_preprocessed_20211020_ID_13785_gg2/NASH_taxonomy_filtered.gg2.asv.counts-rclr.txt"
data_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/"
results_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/stool_mtb/"

#########################################################################
#functions

birdman_dat<-function(dat,annot,suddat){
  df<-dat%>%
    dplyr::rename(ratio=`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_mean`,
                  FeatureID=Feature)%>%
    mutate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`=gsub("[(]|[)]","",`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`))%>%
    separate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`,c("min","max"), sep=",")%>%
    mutate(min=as.numeric(min),
           max=as.numeric(max),
           cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
    left_join(.,annot, by="FeatureID")%>%
    mutate(Name=paste(FeatureID,Compound_Name,mtb_group,sep=" "))%>%
    mutate(Name=ifelse(is.na(Compound_Name),paste(FeatureID, "m/z",mz,mtb_group,sep=" "),Name))
   
  if(suddat=="just_credible"){
    df<-df%>%
      filter(cred=="credible")%>%
      dplyr::select(FeatureID, ratio, min, max, Name)%>%
      arrange(ratio)
  }
  else{
    df<-df%>%
      dplyr::select(FeatureID,Name,mz,mtb_group,ratio, min, max)%>%
      arrange(ratio)
  }
  
  return(df)
}

birdman_dat2<-function(dat,annot,annot2){
  df<-dat%>%
    dplyr::rename(ratio=`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_mean`,
                  FeatureID=Feature)%>%
    mutate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`=gsub("[(]|[)]","",`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`))%>%
    separate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`,c("min","max"), sep=",")%>%
    mutate(min=as.numeric(min),
           max=as.numeric(max),
           cred=ifelse(min>0|max< 0,"credible","not_credible"))%>%
    left_join(.,annot, by="FeatureID")%>%
    mutate(Name=paste(FeatureID,Compound_Name,mtb_group,sep=" "))%>%
    mutate(Name=ifelse(is.na(Compound_Name),paste(FeatureID, "m/z",mz,mtb_group,sep=" "),Name))%>%
    left_join(.,annot2,by="FeatureID")%>%
    dplyr::select(Name,FeatureID,mz,RT,Compound_Name.x,LibraryQualityString,SharedPeaks,LibraryName,
                  MQScore,MassDiff,MZErrorPPM,Adduct,IonMode,mtb_group,ratio, min, max,cred)%>%
    arrange(ratio)
  return(df)
}


birdman_plot<-function(dat){
  dat$Name <- factor(dat$Name,levels = dat$Name)
  
  p<-ggplot(dat, aes(x =Name , y = ratio, ymin = min, ymax = max)) + 
    geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
    geom_pointrange(position = position_dodge(width = 0.8)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip()
  return(p)
}

paired_plots_ind<-function(dat,mtbid){
  count_dat<-dat%>%
    #filter(FeatureID==mtbid & ATTRIBUTE_timepoint=="ZT1" & 
   filter(FeatureID==mtbid & ATTRIBUTE_timepoint=="ZT13" & 
                          ATTRIBUTE_host_subject_id!="19R" & NASH_category!="not_applicable")%>%
    arrange(ATTRIBUTE_disease_stage,ATTRIBUTE_host_subject_id)%>%  
    mutate(log_counts=log10(counts+1))%>%as.data.table()
  
  p<-ggplot(count_dat, aes(x = NASH_category, y = log_counts)) + 
    geom_boxplot(aes(fill = NASH_category), alpha = .5) +
    #geom_line(aes(group = ATTRIBUTE_host_subject_id)) + 
    ggtitle(stringr::str_wrap(unique(count_dat$label_name), width=30))+
    #geom_point(size = 2) + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=2)+
    theme_bw()+ theme(legend.position = "top")+
    scale_fill_manual(values=c("#0000a7","#c1272d"))+ 
    facet_wrap(~ ATTRIBUTE_disease_stage)
  return(p)
}

change_baseline_dat<-function(dat,group){
  df<-dat%>%filter(mtb_group==group)%>%
    mutate(wk128_ratio=log(NASH/(TRF+1)))%>%
    #filter(TRF>0)%>%
    mutate(nash_cond=paste(ATTRIBUTE_condition,NASH_category,sep="_"))
  return(df)
}

change_baseline_plt<-function(dat,type,annot){
  if(type=="combFAFT"){
    p <- ggplot(dat, aes(x=NASH_category, y=wk128_ratio, fill=NASH_category)) + 
      geom_violin(trim=FALSE,alpha=0.7)  + 
      geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), size=0.001)+
      theme_classic() + 
      theme(legend.position = "top") +
      scale_fill_manual(values=c("#0000a7","#c1272d","grey90")) + 
      labs(title=annot,x="NASH Category",y="log2(Wk12/Wk8)")
    
  }
  else{
    p <- ggplot(dat, aes(x=ATTRIBUTE_condition, y=wk128_ratio, fill=NASH_category)) + 
      geom_violin(trim=FALSE,alpha=0.7)  + 
      geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), size=0.01)+
      theme_classic() + 
      theme(legend.position = "top") +
      scale_fill_manual(values=c("#0000a7","#c1272d","grey90")) + 
      labs(title=annot,x="condition",y="log2(Wk12/Wk8)")
  }
  return(p)
}

run_ind_wlcx<-function(data){
    x<-pairwise.wilcox.test(data$wk128_ratio, data$NASH_category, p.adjust.method="fdr")
    df <- data.frame(pval = x$p.value[[1]])
  return(df)
}

run_ind_wlcx_sep<-function(data){
  x<-pairwise.wilcox.test(data$wk128_ratio, data$nash_cond, p.adjust.method="fdr")
  df <- data.frame(pval_FA = x$p.value[[1]], pval_FT=x$p.value[[11]])
  return(df)
}

indv_changebaseline_plt <- function(data,FeatureID) {
  
  lab<-unique((data%>%filter(FeatureID==FeatureID))$label_name)
  
  data%>%
    ggplot(aes(x=NASH_category, y=wk128_ratio, fill=NASH_category)) + 
    geom_violin(alpha=0.7)  + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), size=1.5)+
    theme_classic() + 
    theme(legend.position = "top") +
    scale_fill_manual(values=c("#0000a7","#c1272d","grey90")) + 
    labs(title=stringr::str_wrap(lab, width=30),x="condition",y="log2(Wk12/Wk8)")
  #ggsave(paste("mtb/stool/changebaseline/BA_indv_violins/",FeatureID,".pdf",sep=""), height=4, width=4)
  #ggsave(paste("mtb/stool/changebaseline/AA_indv_violins/",FeatureID,".pdf",sep=""), height=4, width=4)
  #ggsave(paste("mtb/stool/changebaseline/lip_indv_violins/",FeatureID,".pdf",sep=""), height=4, width=4)
  ggsave(paste("mtb/stool/changebaseline/bil_indv_violins/",FeatureID,".pdf",sep=""), height=4, width=4)
}

#########################################################################
#Combined BIRDMAn results--Table S2

annot<-fread("mtb/stool/data_files/20231108_gnps_run/DB_result/8dfc7fa4496841558d688c25f58c63c0.tsv")%>%
  dplyr::rename(FeatureID=`#Scan#`)%>%
  dplyr::select(FeatureID,Precursor_MZ,RT_Query,Compound_Name,LibraryQualityString,SharedPeaks,LibraryName,
                MQScore,MassDiff,MZErrorPPM,Adduct,IonMode)

ZT1<-fread("mtb/stool/birdman/sepZTs/quantification_table-00000-clean-mash-Wk8-ZT1FAFT.beta_var.tsv")
ZT1<-birdman_dat2(ZT1,annotations,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")),
         collection_time="ZT1",disease_stage="Wk 4")

ZT13<-fread("mtb/stool/birdman/sepZTs/quantification_table-00000-clean-mash-Wk8-ZT13FAFT.beta_var.tsv")
ZT13<-birdman_dat2(ZT13,annotations,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")),
         collection_time="ZT13",disease_stage="Wk 4")

ZT12<-fread("mtb/stool/birdman/sepZTs/quantification_table-00000-clean-mash-Wk12-ZT1FAFT.beta_var.tsv")
ZT12<-birdman_dat2(ZT12,annotations,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")),
         collection_time="ZT1",disease_stage="Wk 7")

ZT132<-fread("mtb/stool/birdman/sepZTs/quantification_table-00000-clean-mash-Wk12-ZT13FAFT.beta_var.tsv")
ZT132<-birdman_dat2(ZT132,annotations,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")),
         collection_time="ZT13",disease_stage="Wk 7")

all_res<-rbind(ZT1,ZT13,ZT12,ZT132)
write.table(all_res,"mtb/stool/birdman/sepZTs/birdman_nonNASHNASH_allresults.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
#########################################################################

#birdman
mz_info<-fread('mtb/stool/data_files/20231108_gnps_run/DB_result/mz_annot_v2.txt')

annotations<-fread("mtb/stool/data_files/20231108_gnps_run/DB_result/DB_results_clean.txt")%>%
  select(FeatureID,Compound_Name)%>%
  mutate(FeatureID=as.integer(FeatureID),
         Compound_Name=gsub("\\(predic.*","",Compound_Name),
         Compound_Name=gsub('"',"",Compound_Name))%>%
  right_join(.,mz_info,by="FeatureID")
write.table(annotations,"mtb/stool/data_files/20231108_gnps_run/DB_result/feature_annotations_clean.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#Wk8 ZT1 129 cred
ZT1<-fread("mtb/stool/birdman/sepZTs/quantification_table-00000-clean-mash-Wk8-ZT1FAFT.beta_var.tsv")
ZT1<-birdman_dat(ZT1,annotations)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))%>%
  filter(!grepl("m/z",Name))#%>%
  #filter(FeatureID!=1005)
write.table(ZT1,"mtb/stool/birdman/sepZTs/birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
write.table(ZT1,"mtb/stool/birdman/sepZTs/birdman_justcred_wk8ZT1_nonNASHNASH_results_justannot.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(FT,"mtb/stool/birdman/birdman_justcred_wk8FT_nonNASHNASH_results0.5cred.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT1) + labs(title="ZT1 Wk8\nnon-NASH vs. NASH")
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk8ZT1_nonNASHNASH.pdf",height=20, width=22)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk8ZT1_nonNASHNASH_justannot.pdf",height=6, width=22)

#Wk8 ZT13 204
ZT13<-fread("mtb/stool/birdman/sepZTs/quantification_table-00000-clean-mash-Wk8-ZT13FAFT.beta_var.tsv")
ZT13<-birdman_dat(ZT13,annotations)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))%>%
  filter(!grepl("m/z",Name))
write.table(ZT13,"mtb/stool/birdman/sepZTs/birdman_justcred_wk8ZT13_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
write.table(ZT13,"mtb/stool/birdman/sepZTs/birdman_justcred_wk8ZT13_nonNASHNASH_results_justannot.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT13) + labs(title="ZT13 Wk8\nnon-NASH vs. NASH")
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk8ZT13_nonNASHNASH.pdf",height=26, width=22)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk8ZT13_nonNASHNASH_justannot.pdf",height=8, width=22)

#Wk12 ZT1 466
ZT1<-fread("mtb/stool/birdman/sepZTs/quantification_table-00000-clean-mash-Wk12-ZT1FAFT.beta_var.tsv")
ZT1<-birdman_dat(ZT1,annotations)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))%>%
  filter(!grepl("m/z",Name))
write.table(ZT1,"mtb/stool/birdman/sepZTs/birdman_justcred_wk12ZT1_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
write.table(ZT1,"mtb/stool/birdman/sepZTs/birdman_justcred_wk12ZT1_nonNASHNASH_results_justannot.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(FA,"mtb/stool/birdman/birdman_justcred_wk8FA_nonNASHNASH_results0.5cred.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT1) + labs(title="ZT1 Wk12\nnon-NASH vs. NASH")
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk12ZT1_nonNASHNASH.pdf",height=25, width=21)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk12ZT1_nonNASHNASH_justannot.pdf",height=12, width=22)

#Wk12 ZT13 634
ZT13<-fread("mtb/stool/birdman/sepZTs/quantification_table-00000-clean-mash-Wk12-ZT13FAFT.beta_var.tsv")
ZT13<-birdman_dat(ZT13,annotations)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))%>%
  filter(!grepl("m/z",Name))
write.table(ZT13,"mtb/stool/birdman/sepZTs/birdman_justcred_wk12ZT13_nonNASHNASH_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
write.table(ZT13,"mtb/stool/birdman/sepZTs/birdman_justcred_wk12ZT13_nonNASHNASH_results_justannot.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

birdman_plot(ZT13) + labs(title="ZT13 Wk12\nnon-NASH vs. NASH")
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk12ZT13_nonNASHNASH.pdf",height=25, width=21)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk12ZT13_nonNASHNASH_justannot.pdf",height=18, width=20)

#########################################################################
##ZT1 Wk8 to Wk12

wk8ZT1<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv")
wk12ZT1<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk12ZT1_nonNASHNASH_results.tsv")

list_venn <- list(Wk8 = wk8ZT1$Name,
                  Wk12 = wk12ZT1$Name)

p<-venn.diagram(list_venn, 
                main="ZT1 Non-NASH vs. NASH",filename = NULL)

grid.draw(p)
pdf(file="mtb/stool/birdman/sepZTs/birdman_justcredoverlap_nonNASHNASH_ZT1Wk8Wk12results.pdf")
grid.draw(p)
dev.off()


##ZT13 Wk8 to Wk12

wk8ZT13<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv")
wk12ZT13<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk12ZT13_nonNASHNASH_results.tsv")

list_venn <- list(Wk8 = wk8ZT13$Name,
                  Wk12 = wk12ZT13$Name)

p<-venn.diagram(list_venn, 
                main="ZT13 non-NASH vs. NASH",filename = NULL)

grid.draw(p)
pdf(file="mtb/stool/birdman/sepZTs/birdman_justcredoverlap_nonNASHNASH_ZT13Wk8Wk12results.pdf")
grid.draw(p)
dev.off()


#comparing all Wk8 and wk12 to ZT1 and ZT13
list_venn <- list(Wk8_Z1 = wk8ZT1$Name,
                  Wk12_Z1 = wk12ZT1$Name,
                  Wk8_Z13 = wk8ZT13$Name,
                  Wk12_Z13 = wk12ZT13$Name)

p<-venn.diagram(list_venn,height = 10,
                main="NASH",
                width = 10,alpha = c(0.5, 0.5, 0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="mtb/stool/birdman/sepZTs/birdman_justcredoverlap_nonNASHNASH_wk8wk12_results.pdf")
grid.draw(p)
dev.off()

##Wk8 ZT1 vs ZT13

wk8ZT1<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv")
wk8ZT13<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk8ZT13_nonNASHNASH_results.tsv")

list_venn <- list(ZT1 = wk8ZT1$Name,
                  ZT13 = wk8ZT13$Name)

p<-venn.diagram(list_venn, fill = c("snow2","gray25"),height = 10,
                main="NASH Wk8",
                width = 10,alpha = c(0.5, 0.5),filename = NULL)

grid.draw(p)
pdf(file="mtb/stool/birdman/sepZTs/birdman_justcredoverlap_nonNASHNASH_wk8ZT113results.pdf")
grid.draw(p)
dev.off()


##Wk12 FA vs FT 

wk12ZT1<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk12ZT1_nonNASHNASH_results.tsv")
wk12ZT13<-fread("mtb/stool/birdman/sepZTs/birdman_justcred_wk12ZT13_nonNASHNASH_results.tsv")

list_venn <- list(ZT1 = wk12ZT1$Name,
                  ZT13 = wk12ZT13$Name)

p<-venn.diagram(list_venn, fill = c("snow2","gray25"),height = 10,
                main="NASH Wk12",
                width = 10,alpha = c(0.5, 0.5), filename = NULL)

grid.draw(p)
pdf(file="mtb/stool/birdman/sepZTs/birdman_justcredoverlap_nonNASHNASH_wk12ZT113results.pdf")
grid.draw(p)
dev.off()

#########################################################################

#plot just the mtb sig hits that overlapped b/w wk8 and wk12

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
intersect(all1$`Wk8:Wk12`,all2$`Wk8:Wk12`) #only 3 intersect 

mtb_list<-c(all1$`Wk8:Wk12`,all2$`Wk8:Wk12`)%>%unique() #41 unique 


#wk8
ZT1<-wk8ZT1%>%filter(Name %in% all1$`Wk8:Wk12`)%>%filter(ratio>1 | ratio< -1)
birdman_plot(ZT1) + labs(title="ZT1 Wk8\nnon-NASH vs. NASH")
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk8ZT1_shared29_nonNASHNASH.pdf",height=6, width=22)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk8ZT1_shared24_1rat_nonNASHNASH.pdf",height=5, width=7)

ZT13<-wk8ZT13%>%filter(Name %in% all2$`Wk8:Wk12`)%>%filter(ratio>1 | ratio< -1)
birdman_plot(ZT13) + labs(title="ZT13 Wk8\nnon-NASH vs. NASH")
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk8ZT13_shared15_nonNASHNASH.pdf",height=4, width=21)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk8ZT13_shared10_1rat_nonNASHNASH.pdf",height=3 , width=21)

#wk12
ZT1<-wk12ZT1%>%filter(Name %in% all1$`Wk8:Wk12`)%>%filter(ratio>1 | ratio< -1)
birdman_plot(ZT1) + labs(title="ZT1 Wk12\nnon-NASH vs. NASH")
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk12ZT1_shared29_nonNASHNASH.pdf",height=6, width=22)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk12ZT1_shared23_1rat_nonNASHNASH.pdf",height=5, width=10)

ZT13<-wk12ZT13%>%filter(Name %in% all2$`Wk8:Wk12`)%>%filter(ratio>1 | ratio< -1)
birdman_plot(ZT13) + labs(title="ZT13 Wk12\nnon-NASH vs. NASH")
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk12ZT13_shared15_nonNASHNASH.pdf",height=4, width=21)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_Wk12ZT13_shared14_1rat_nonNASHNASH.pdf",height=4 , width=21)

#ZT1
w8<-wk8ZT1%>%filter(Name %in% all1$`Wk8:Wk12`)%>%
  mutate(grp="WK8")%>%arrange(ratio)
w12<-wk12ZT1%>%filter(Name %in% all1$`Wk8:Wk12`)%>%
  mutate(grp="Wk12")%>%arrange(ratio)

combZT1<-rbind(w8,w12)

combZT1$Name <- factor(combZT1$Name,levels = w12$Name)
#combZT1$Name <- factor(combZT1$Name,levels = w8$Name)

p<-ggplot(combZT1, aes(x =Name , y = ratio, ymin = min, ymax = max, color=grp)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_color_manual(values=c("black","gray60")) + ggtitle("ZT1")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip()
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_shared29_nonNASHNASH.pdf",height=6 , width=24)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_shared29_ordwk8_nonNASHNASH.pdf",height=6 , width=24)

#ZT13
w8<-wk8ZT13%>%filter(Name %in% all2$`Wk8:Wk12`)%>%
  mutate(grp="WK8")%>%arrange(ratio)
w12<-wk12ZT13%>%filter(Name %in% all2$`Wk8:Wk12`)%>%
  mutate(grp="Wk12")%>%arrange(ratio)

combZT13<-rbind(w8,w12)

combZT13$Name <- factor(combZT13$Name,levels = w12$Name)
#combZT13$Name <- factor(combZT13$Name,levels = w8$Name)

p<-ggplot(combZT13, aes(x =Name , y = ratio, ymin = min, ymax = max, color=grp)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_color_manual(values=c("black","gray60")) + ggtitle("ZT13")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip()
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT13_wk812comb_shared15_nonNASHNASH.pdf",height=3.5 , width=23)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT13_wk812comb_shared15_ordwk8_nonNASHNASH.pdf",height=3.5 , width=22)

#plot specific examples

mtb<-fread("mtb/stool/NASH_stool_metabolomics_quanttab_long_wmd_v2.txt")%>%
  dplyr::select(c(1,4:5,3,11,21:23,27,33))%>%
  mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition),
         cln_Compound_Name=gsub("\\(predicted.*","",Compound_Name))%>%
  mutate(cln_Compound_Name=gsub("\\(delta.*","",cln_Compound_Name))%>%
  mutate(label_name=ifelse(is.na(cln_Compound_Name),paste("mtb", FeatureID, "m/z",mz, sep=" "),paste("mtb",FeatureID,"m/z",mz, cln_Compound_Name, sep=" ")))%>%
  mutate(ATTRIBUTE_disease_stage=factor(ATTRIBUTE_disease_stage,levels=c("TRF","NASH")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","not_applicable")),
         ATTRIBUTE_condition=factor(ATTRIBUTE_condition,levels=c("NA","FA","FT")))
mtb1264<-paired_plots_ind(mtb,1264)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_mtb1264_nonNASHNASH.pdf",height=4 , width=4)

mtb3597<-paired_plots_ind(mtb,3597)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_mtb3597_nonNASHNASH.pdf",height=4 , width=4)

mtb1992<-paired_plots_ind(mtb,1992)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_mtb1992_nonNASHNASH.pdf",height=4 , width=4)

mtb3748<-paired_plots_ind(mtb,3748)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_mtb3748_nonNASHNASH.pdf",height=4 , width=4)

mtb3552<-paired_plots_ind(mtb,3552)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_mtb3552_nonNASHNASH.pdf",height=4 , width=4)

mtb2245<-paired_plots_ind(mtb,2245)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_mtb2245_nonNASHNASH.pdf",height=4 , width=4)

mtb2568<-paired_plots_ind(mtb,2568)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_mtb2568_nonNASHNASH.pdf",height=4 , width=4)

mtb2691<-paired_plots_ind(mtb,2691)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_mtb2691_nonNASHNASH.pdf",height=4 , width=4)

mtb4446<-paired_plots_ind(mtb,4446)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_mtb4446_nonNASHNASH.pdf",height=4 , width=4)

##biological relevant ones
mtb1453<-paired_plots_ind(mtb,1453)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_mtb1453_nonNASHNASH.pdf",height=4 , width=4)

mtb2908<-paired_plots_ind(mtb,2908)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT1_wk812comb_mtb2908_nonNASHNASH.pdf",height=4 , width=4)

#ZT13
mtb2245<-paired_plots_ind(mtb,2245)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT13_wk812comb_mtb2245_nonNASHNASH.pdf",height=4 , width=4)

mtb4168<-paired_plots_ind(mtb,4168)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT13_wk812comb_mtb4168_nonNASHNASH.pdf",height=4 , width=4)

mtb4200<-paired_plots_ind(mtb,4200)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT13_wk812comb_mtb4200_nonNASHNASH.pdf",height=4 , width=4)

mtb4166<-paired_plots_ind(mtb,4166)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT13_wk812comb_mtb4166_nonNASHNASH.pdf",height=4 , width=4)

mtb1522<-paired_plots_ind(mtb,1522)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT13_wk812comb_mtb1522_nonNASHNASH.pdf",height=4 , width=4)

##biological relevant ones
mtb1453<-paired_plots_ind(mtb,1453)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT13_wk812comb_mtb1453_nonNASHNASH.pdf",height=4 , width=4)

mtb2908<-paired_plots_ind(mtb,2908)
ggsave("mtb/stool/birdman/sepZTs/birdman_justcred_ZT13_wk812comb_mtb2908_nonNASHNASH.pdf",height=4 , width=4)

#########################################################################

#take the log of wk12/wk8 for BA non-nash vs. nash

# class<-fread("mtb/stool/old_files/canopus_summary.tsv")%>%
#   separate(name, c(NA,"FeatureID"),sep="sirius_")%>%
#   mutate(FeatureID=as.integer(FeatureID))
#feat<-fread("mtb/stool/old_files/feature_annotations_clean.txt")
# md<-fread("mtb/stool/old_files/NASH_stool_metabolomics_metadata.txt")%>%
#   mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition))%>%
#   separate(sample_id_mouse, c("mouse_id",NA))

#annot<-fread("mtb/stool/data_files/20231108_gnps_run/DB_result/DB_results_clean.txt")
#mz<-fread("mtb/stool/data_files/20231108_gnps_run/DB_result/mz_annot.txt")

#check for mice that had both wk8 and wk12 samples
wk8wk12samps<-fread("mtb/stool/NASH_stool_metabolomics_metadata_new.txt")%>%
  mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition))%>%
  filter(!is.na(NASH_category))%>%
  group_by(ATTRIBUTE_host_subject_id,ATTRIBUTE_condition,ATTRIBUTE_timepoint)%>%
  summarise(n=n())%>%
  arrange(ATTRIBUTE_host_subject_id)%>%
  filter(n>1)

# mtb<-read_qza("files_from_KF/STAM-TRF-ileum-cecum-16S-study-13785/metabolomics/stool/data/filtered_featuretable_allzt_alldiets_NASH.qza")$data%>%
#   as.data.frame()%>%rownames_to_column("FeatureID")%>%
#mtb<-fread("mtb/stool/NASH_stool_metabolomics_quanttab_long_wmd.txt")%>%
mtb<-fread("mtb/stool/NASH_stool_metabolomics_quanttab_long_wmd_v2.txt")%>%
  dplyr::select(c(1,4:5,3,11,21:23,27,33))%>%
  mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition))%>%
  filter(ATTRIBUTE_host_subject_id %in% wk8wk12samps$ATTRIBUTE_host_subject_id)%>%
  filter(ATTRIBUTE_timepoint=="ZT1")%>%
  spread(ATTRIBUTE_disease_stage,counts)%>%
  arrange(ATTRIBUTE_host_subject_id)%>%
  filter(!is.na(NASH))%>%
  filter(!is.na(TRF))%>%
  mutate(NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","not_applicable")),
         ATTRIBUTE_condition=factor(ATTRIBUTE_condition,levels=c("NA","FA","FT")))

mtb_dens<-mtb%>%
  mutate(wk128_ratio=log(NASH/(TRF+1)))%>%
  mutate(nash_cond=paste(ATTRIBUTE_condition,NASH_category,sep="_"))%>%
  filter(mtb_group!="unknown")%>% #& ATTRIBUTE_condition!="NA"
 filter(mtb_group!="glycols")

p <- ggplot(mtb_dens, aes(x=wk128_ratio, colour=NASH_category)) + theme_classic()+
  geom_density() + scale_colour_manual(values=c("#0000a7","#c1272d","grey70"))+ 
  geom_vline(aes(xintercept=0),linetype="dashed") +
  facet_grid(ATTRIBUTE_condition~mtb_group)

ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_density_wNA.pdf",p,height=4, width=12)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_woglycols_density_wNA.pdf",p,height=4, width=10)

p <- ggplot(mtb_dens, aes(x=wk128_ratio, colour=NASH_category)) + theme_classic()+
  geom_density() + scale_colour_manual(values=c("#0000a7","#c1272d","grey70"))+ 
  geom_vline(aes(xintercept=0),linetype="dashed") +
  facet_grid(~mtb_group)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_woglycols_density_wNA_combFAFT.pdf",p,height=2, width=10)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_density_wNA_combFAFT.pdf",p,height=2, width=12)

#bile acids
BA_mtb<-change_baseline_dat(mtb,"bile acids")%>%
  filter(ATTRIBUTE_condition!="NA")
write.table(BA_mtb,"mtb/stool/changebaseline/foldchange_wk12wk8_BAs.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

p <- change_baseline_plt(BA_mtb,"normal","Bile Acids")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_BA.pdf",p,height=4, width=5)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_BA_woNA.pdf",p,height=4, width=5)

pairwise.wilcox.test(BA_mtb$wk128_ratio, BA_mtb$nash_cond,
                     p.adjust.method="fdr")

#                 FA_NASH FA_Non_NASH FT_NASH FT_Non_NASH
# FA_Non_NASH       0.980   -           -       -          
#   FT_NASH           3.7e-13 5.7e-14     -       -          
#   FT_Non_NASH       1.3e-14 6.3e-15     0.047   -          
#   NA_not_applicable 1.0e-14 6.4e-15     0.115   0.558   

p <- change_baseline_plt(BA_mtb,"combFAFT","Bile Acids")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_BA_combFAFT.pdf",p,height=4, width=4)

pairwise.wilcox.test(BA_mtb$wk128_ratio, BA_mtb$NASH_category,
                     p.adjust.method="fdr")
#             Non_NASH NASH   
# NASH           0.66     -      
#   not_applicable 5.4e-06  2.0e-06

#run wilcox on indiv mtb to find non-nash vs nash diff
BA_inds<-BA_mtb%>%
  mutate(cln_Compound_Name=gsub("\\(predicted.*","",Compound_Name))%>%
  mutate(label_name=ifelse(is.na(cln_Compound_Name),paste("mtb", FeatureID, "m/z",mz, sep=" "),paste("mtb",FeatureID,"m/z",mz, cln_Compound_Name, sep=" ")))

BA_nested <- BA_inds %>% 
  group_by(FeatureID) %>% 
  nest()

mtbBA_pval <- 
  BA_nested %>% 
  mutate(pval = map2(data,FeatureID,  ~ run_ind_wlcx(.x)),
         pval_sep=map2(data,FeatureID,  ~ run_ind_wlcx_sep(.x)))%>%
  dplyr::select(pval,pval_sep)%>%
  unnest() #3 were borderline different

write.table(mtbBA_pval,"mtb/stool/changebaseline/foldchange_wk12wk8_BAs_pvals.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

#individual BA plots (not separating by FA or FT due to power)
mtbBA_plots <- 
  BA_nested %>% 
  mutate(plot = map2(data, FeatureID,  ~ indv_changebaseline_plt(.x,.y)))

order<-BA_mtb%>%
  mutate(FeatureID=as.factor(FeatureID))%>%
  #filter(ATTRIBUTE_condition=="FT" & wk128_ratio!=-Inf & NASH_category=="Non_NASH")%>%
  #filter(ATTRIBUTE_condition=="FA" & wk128_ratio!=-Inf & NASH_category=="Non_NASH")%>%
  filter(ATTRIBUTE_condition=="NA" & wk128_ratio!=-Inf)%>%
  group_by(FeatureID,ATTRIBUTE_condition,NASH_category)%>%
  summarise(mean_ratio=mean(wk128_ratio))%>%
  ungroup()%>%
  arrange(mean_ratio)%>%
  mutate(FeatureID=fct_reorder(FeatureID, mean_ratio))

BAranked<-BA_mtb %>%
  filter(wk128_ratio!=-Inf) %>%
  group_by(FeatureID,ATTRIBUTE_condition,NASH_category) %>%
  summarise(mean_ratio=mean(wk128_ratio)) %>%
  ungroup() %>%
  mutate(FeatureID = factor(FeatureID,levels=order$FeatureID))


p1<-ggplot(data=BAranked, aes(x=FeatureID, y=mean_ratio, fill=NASH_category)) +
  geom_bar(stat="identity", position=position_dodge()) +scale_fill_manual(values=c("#0000a7","#c1272d","grey70")) +
  facet_grid(rows = vars(ATTRIBUTE_condition)) +theme_classic()+
  labs(title="Bile acids", x="Features",y="mean log2(Wk12/Wk8)")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_BA_ranked_byNA.pdf",p1,height=4, width=10)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_BA_ranked_byFAnonnash.pdf",p1,height=4, width=10)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_BA_ranked_byFTnonnash.pdf",p1,height=4, width=10)

#combined rank plot
order<-BA_mtb%>%
  mutate(FeatureID=as.factor(FeatureID))%>%
  #filter(wk128_ratio!=-Inf & NASH_category=="Non_NASH")%>%
  filter(ATTRIBUTE_condition=="NA" & wk128_ratio!=-Inf)%>%
  group_by(FeatureID,NASH_category)%>%
  summarise(mean_ratio=mean(wk128_ratio))%>%
  ungroup()%>%
  arrange(mean_ratio)%>%
  mutate(FeatureID=fct_reorder(FeatureID, mean_ratio))

BAranked<-BA_mtb %>%
  filter(wk128_ratio!=-Inf) %>%
  group_by(FeatureID,NASH_category) %>%
  summarise(mean_ratio=mean(wk128_ratio)) %>%
  ungroup() %>%
  mutate(FeatureID = factor(FeatureID,levels=order$FeatureID))

p1<-ggplot(data=BAranked, aes(x=FeatureID, y=mean_ratio, fill=NASH_category)) +
  geom_bar(stat="identity", position=position_dodge()) +scale_fill_manual(values=c("#0000a7","#c1272d","grey70")) +
  theme_classic()+
  labs(title="Bile acids", x="Features",y="mean log2(Wk12/Wk8)")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_BA_ranked_combbyNA.pdf",p1,height=2.5, width=10)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_BA_ranked_combbynonnash.pdf",p1,height=2.5, width=10)

#amino acids
AA_mtb<-change_baseline_dat(mtb,"amines")%>%
  filter(ATTRIBUTE_condition!="NA")
write.table(AA_mtb,"mtb/stool/changebaseline/foldchange_wk12wk8_AAs.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

p <- change_baseline_plt(AA_mtb,"normal","Amino acids, peptides, and analogues")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_AA.pdf",p,height=4, width=5)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_AA_woNA.pdf",p,height=4, width=5)

pairwise.wilcox.test(AA_mtb$wk128_ratio, AA_mtb$nash_cond,
                     p.adjust.method="fdr")

#                  FA_NASH FA_Non_NASH FT_NASH FT_Non_NASH
# FA_Non_NASH       <2e-16  -           -       -          
#   FT_NASH           0.035   <2e-16      -       -          
#   FT_Non_NASH       0.597   <2e-16      0.096   -          
#   NA_not_applicable <2e-16  <2e-16      <2e-16  <2e-16 

p <- change_baseline_plt(AA_mtb,"combFAFT","Amino acids, peptides, and analogues")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_AA_combFAFT.pdf",p,height=4, width=4)

pairwise.wilcox.test(AA_mtb$wk128_ratio, AA_mtb$NASH_category,
                     p.adjust.method="fdr")

#               Non_NASH NASH  
# NASH           <2e-16   -     
#   not_applicable <2e-16   <2e-16

#sun wilcox on indiv mtb to find non-nash vs nash diff
AA_inds<-AA_mtb%>%
  mutate(cln_Compound_Name=gsub("\\(predicted.*","",Compound_Name))%>%
  mutate(label_name=ifelse(is.na(cln_Compound_Name),paste("mtb", FeatureID, "m/z",mz, sep=" "),paste("mtb",FeatureID,"m/z",mz, cln_Compound_Name, sep=" ")))

AA_nested <- AA_inds %>% 
  group_by(FeatureID) %>% 
  nest()

mtbAA_pval <- 
  AA_nested %>% 
  mutate(pval = map2(data,FeatureID,  ~ run_ind_wlcx(.x)),
         pval_sep=map2(data,FeatureID,  ~ run_ind_wlcx_sep(.x)))%>%
  dplyr::select(pval,pval_sep)%>%
  unnest() #3 were borderline different

write.table(mtbAA_pval,"mtb/stool/changebaseline/foldchange_wk12wk8_AAs_pvals.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

#individual BA plots (not separating by FA or FT due to power)
mtbAA_plots <- 
  AA_nested %>% 
  mutate(plot = map2(data, FeatureID,  ~ indv_changebaseline_plt(.x,.y)))

#lipids/fatty acids
lip_mtb<-change_baseline_dat(mtb,"fatty acids & other lipids")%>%
  filter(ATTRIBUTE_condition!="NA")
write.table(lip_mtb,"mtb/stool/changebaseline/foldchange_wk12wk8_lip.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

p <- change_baseline_plt(lip_mtb,"normal","fatty acids & other lipids")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_lipid.pdf",p,height=4, width=5)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_lipid_woNA.pdf",p,height=4, width=5)

pairwise.wilcox.test(lip_mtb$wk128_ratio, lip_mtb$nash_cond,
                     p.adjust.method="fdr")

#                   FA_NASH FA_Non_NASH FT_NASH FT_Non_NASH
# FA_Non_NASH       1.6e-10 -           -       -          
#   FT_NASH           0.143   < 2e-16     -       -          
#   FT_Non_NASH       7.8e-06 0.013       6.1e-13 -          
#   NA_not_applicable < 2e-16 < 2e-16     < 2e-16 < 2e-16   

p <- change_baseline_plt(lip_mtb,"combFAFT","fatty acids & other lipids")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_lipid_combFAFT.pdf",p,height=4, width=4)

pairwise.wilcox.test(lip_mtb$wk128_ratio, lip_mtb$NASH_category,
                     p.adjust.method="fdr")

# Non_NASH NASH  
# NASH           <2e-16   -     
#   not_applicable <2e-16   <2e-16

#sun wilcox on indiv mtb to find non-nash vs nash diff
lip_inds<-lip_mtb%>%
  mutate(cln_Compound_Name=gsub("\\(predicted.*","",Compound_Name))%>%
  mutate(label_name=ifelse(is.na(cln_Compound_Name),paste("mtb", FeatureID, "m/z",mz, sep=" "),paste("mtb",FeatureID,"m/z",mz, cln_Compound_Name, sep=" ")))

lip_nested <- lip_inds %>% 
  group_by(FeatureID) %>% 
  nest()

mtblip_pval <- 
  lip_nested %>% 
  mutate(pval = map2(data,FeatureID,  ~ run_ind_wlcx(.x)),
         pval_sep=map2(data,FeatureID,  ~ run_ind_wlcx_sep(.x)))%>%
  dplyr::select(pval,pval_sep)%>%
  unnest() #11 were borderline different

write.table(mtblip_pval,"mtb/stool/changebaseline/foldchange_wk12wk8_lipid_pvals.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

#individual plots (not separating by FA or FT due to power)
mtblip_plots <- 
  lip_nested %>% 
  mutate(plot = map2(data, FeatureID,  ~ indv_changebaseline_plt(.x,.y)))

#glycols
gly_mtb<-change_baseline_dat(mtb,"glycols")%>%
  filter(ATTRIBUTE_condition!="NA")
write.table(gly_mtb,"mtb/stool/changebaseline/foldchange_wk12wk8_gly.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

p <- change_baseline_plt(gly_mtb,"normal","glycols")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_glycol.pdf",p,height=4, width=5)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_glycol_woNA.pdf",p,height=4, width=5)

pairwise.wilcox.test(gly_mtb$wk128_ratio, gly_mtb$nash_cond,
                     p.adjust.method="fdr")

# FA_NASH FA_Non_NASH FT_NASH FT_Non_NASH
# FA_Non_NASH       0.05608 -           -       -          
#   FT_NASH           0.51433 0.00498     -       -          
#   FT_Non_NASH       0.09949 0.00069     0.13222 -          
#   NA_not_applicable 0.91534 0.01415     0.41594 0.03651   

p <- change_baseline_plt(gly_mtb,"combFAFT","glycols")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_glycol_combFAFT.pdf",p,height=4, width=4)

pairwise.wilcox.test(gly_mtb$wk128_ratio, gly_mtb$NASH_category,
                     p.adjust.method="fdr")
# Non_NASH NASH
# NASH           0.75     -   
#   not_applicable 0.79     0.75

#carnitines
car_mtb<-change_baseline_dat(mtb,"carnitines")%>%
  filter(ATTRIBUTE_condition!="NA")
write.table(car_mtb,"mtb/stool/changebaseline/foldchange_wk12wk8_car.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

p <- change_baseline_plt(car_mtb,"normal","carnitines")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_carnitines.pdf",p,height=4, width=5)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_carnitines_woNA.pdf",p,height=4, width=5)

pairwise.wilcox.test(car_mtb$wk128_ratio, car_mtb$nash_cond,
                     p.adjust.method="fdr")

# FA_NASH FA_Non_NASH FT_NASH FT_Non_NASH
# FA_Non_NASH       0.92    -           -       -          
#   FT_NASH           0.83    0.83        -       -          
#   FT_Non_NASH       0.91    0.91        0.59    -          
#   NA_not_applicable 0.92    0.91        0.61    0.91  

p <- change_baseline_plt(car_mtb,"combFAFT","carnitines")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_carnitines_combFAFT.pdf",p,height=4, width=4)

pairwise.wilcox.test(car_mtb$wk128_ratio, car_mtb$NASH_category,
                     p.adjust.method="fdr")

# Non_NASH NASH
# NASH           0.34     -   
#   not_applicable 0.94     0.34

#phosphocholines
phos_mtb<-change_baseline_dat(mtb,"phosphocholines")%>%
  filter(ATTRIBUTE_condition!="NA")
write.table(phos_mtb,"mtb/stool/changebaseline/foldchange_wk12wk8_phos.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

p <- change_baseline_plt(phos_mtb,"normal","phosphocholines")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_phosphocholine.pdf",p,height=4, width=5)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_phosphocholine_woNA.pdf",p,height=4, width=5)

pairwise.wilcox.test(phos_mtb$wk128_ratio, phos_mtb$nash_cond,
                     p.adjust.method="fdr")

#                  FA_NASH FA_Non_NASH FT_NASH FT_Non_NASH
# FA_Non_NASH       0.0024  -           -       -          
#   FT_NASH           0.1139  0.0090      -       -          
#   FT_Non_NASH       0.0268  6.2e-10     8.6e-08 -          
#   NA_not_applicable 0.2691  8.1e-05     0.0202  0.0559  

p <- change_baseline_plt(phos_mtb,"combFAFT","phosphocholines")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_phosphocholine_combFAFT.pdf",p,height=4, width=4)

pairwise.wilcox.test(phos_mtb$wk128_ratio, phos_mtb$NASH_category,
                     p.adjust.method="fdr")
             
#             Non_NASH NASH
# NASH           0.31     -   
#   not_applicable 0.24     0.05

#bilirubin
bil_mtb<-change_baseline_dat(mtb,"bilirubins")%>%
  filter(ATTRIBUTE_condition!="NA")
write.table(bil_mtb,"mtb/stool/changebaseline/foldchange_wk12wk8_bil.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

p <- change_baseline_plt(bil_mtb,"normal","bilirubins")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_bilirubins.pdf",p,height=4, width=5)
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_bilirubins_woNA.pdf",p,height=4, width=5)

pairwise.wilcox.test(bil_mtb$wk128_ratio, bil_mtb$nash_cond,
                     p.adjust.method="fdr")

# FA_NASH FA_Non_NASH FT_NASH FT_Non_NASH
# FA_Non_NASH       0.0969  -           -       -          
#   FT_NASH           0.0027  1.6e-09     -       -          
#   FT_Non_NASH       0.1020  0.9886      1.8e-09 -          
#   NA_not_applicable 2.5e-06 8.2e-14     0.0449  8.2e-14   

p <- change_baseline_plt(bil_mtb,"combFAFT","bilirubins")
ggsave("mtb/stool/changebaseline/foldchange_wk12wk8_bilirubin_combFAFT.pdf",p,height=4, width=4)

pairwise.wilcox.test(bil_mtb$wk128_ratio, bil_mtb$NASH_category,
                     p.adjust.method="fdr")
# Non_NASH NASH 
# NASH           1.8e-10  -    
#   not_applicable < 2e-16  1e-04

#sun wilcox on indiv mtb to find non-nash vs nash diff
bil_inds<-bil_mtb%>%
  mutate(cln_Compound_Name=gsub("\\(predicted.*","",Compound_Name))%>%
  mutate(label_name=ifelse(is.na(cln_Compound_Name),paste("mtb", FeatureID, "m/z",mz, sep=" "),paste("mtb",FeatureID,"m/z",mz, cln_Compound_Name, sep=" ")))

bil_nested <- bil_inds %>% 
  group_by(FeatureID) %>% 
  nest()

mtbbil_pval <- 
  bil_nested %>% 
  mutate(pval = map2(data,FeatureID,  ~ run_ind_wlcx(.x)),
         pval_sep=map2(data,FeatureID,  ~ run_ind_wlcx_sep(.x)))%>%
  dplyr::select(pval,pval_sep)%>%
  unnest() #2 were borderline different

write.table(mtbbil_pval,"mtb/stool/changebaseline/foldchange_wk12wk8_bilirubin_pvals.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

#individual plots (not separating by FA or FT due to power)
mtbbil_plots <- 
  bil_nested %>% 
  mutate(plot = map2(data, FeatureID,  ~ indv_changebaseline_plt(.x,.y)))
