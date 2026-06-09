setwd("~/Notebooks/sfloresr/MASH-TRF/MASHomics/")

library(tidyverse)
library(data.table)
library(VennDiagram)
library(gplots)

##################################################################
#inputs
mash_metadata<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/metadata_cln_addmoreNASHcat.txt"
mtb_annotations<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/20231108_gnps_run/DB_result/8dfc7fa4496841558d688c25f58c63c0.tsv"
mtb_clean<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/20231108_gnps_run/DB_result/feature_annotations_clean.txt"
bdm_res_w8_ZT1<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/stool_mtb/birdman_results/quantification_table-00000-clean-mash-Wk8-ZT1FAFT.beta_var.tsv"
bdm_res_w8_ZT13<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/stool_mtb/birdman_results/quantification_table-00000-clean-mash-Wk8-ZT13FAFT.beta_var.tsv"
bdm_res_w12_ZT1<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/stool_mtb/birdman_results/quantification_table-00000-clean-mash-Wk12-ZT1FAFT.beta_var.tsv"
bdm_res_w12_ZT13<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/stool_mtb/birdman_results/quantification_table-00000-clean-mash-Wk12-ZT13FAFT.beta_var.tsv"
feattab<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/NASH_stool_metabolomics_quanttab_long_wmd_v2.txt"
data_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_mtb/"
bdm_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/stool_mtb/birdman_results/"

#########################################################################
#functions

birdman_dat<-function(dat,annot,justannot){
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
    filter(cred=="credible") %>%
    dplyr::select(FeatureID,Name,mz,mtb_group,ratio, min, max)%>%
    arrange(ratio)%>%
    mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
    mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
  
  if(justannot=="just_annot"){
    df<-df%>%filter(!grepl("m/z",Name))
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
    geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+#theme_pubr()+
    geom_pointrange(position = position_dodge(width = 0.8)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip()
  return(p)
}

paired_plots_ind<-function(dat,mtbid,tmpt){
  count_dat<-dat
  
  if(tmpt=="ZT1"){
    count_dat<-count_dat%>%
      filter(FeatureID==mtbid & ATTRIBUTE_timepoint=="ZT1" & 
               ATTRIBUTE_host_subject_id!="19R" & NASH_category!="not_applicable")
  }else{
    count_dat<-count_dat%>%
      filter(FeatureID==mtbid & ATTRIBUTE_timepoint=="ZT13" & 
               ATTRIBUTE_host_subject_id!="19R" & NASH_category!="not_applicable")
    }
  
  count_dat<-count_dat%>%
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

#########################################################################
#Combined BIRDMAn results--Table S2

annot<-fread(mtb_annotations)%>%
  dplyr::rename(FeatureID=`#Scan#`)%>%
  dplyr::select(FeatureID,Precursor_MZ,RT_Query,Compound_Name,LibraryQualityString,SharedPeaks,LibraryName,
                MQScore,MassDiff,MZErrorPPM,Adduct,IonMode)

ZT1<-fread(bdm_res_w8_ZT1)
ZT1<-birdman_dat2(ZT1,annotations,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")),
         collection_time="ZT1",disease_stage="Wk 4")

ZT13<-fread(bdm_res_w8_ZT13)
ZT13<-birdman_dat2(ZT13,annotations,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")),
         collection_time="ZT13",disease_stage="Wk 4")

ZT12<-fread(bdm_res_w12_ZT1)
ZT12<-birdman_dat2(ZT12,annotations,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")),
         collection_time="ZT1",disease_stage="Wk 7")

ZT132<-fread(bdm_res_w12_ZT13)
ZT132<-birdman_dat2(ZT132,annotations,annot)%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")),
         collection_time="ZT13",disease_stage="Wk 7")

all_res<-rbind(ZT1,ZT13,ZT12,ZT132)
write.table(all_res,paste0(bdm_dir,"birdman_nonNASHNASH_allresults.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)
#########################################################################
#Differential abundance analysis (BIRDMAn) results NASH vs. non-NASH (FA)

annotations<-fread(mtb_clean)

#Wk8 ZT1 129
ZT1<-fread(bdm_res_w8_ZT1)
ZT1<-birdman_dat(ZT1,annotations,"all")
write.table(ZT1,paste0(bdm_dir,"birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 
birdman_plot(ZT1) + labs(title="ZT1 Wk8\nnon-NASH vs. NASH")
ggsave(paste0(bdm_dir,"birdman_justcred_Wk8ZT1_nonNASHNASH.pdf"),height=20, width=22)

ZT1<-fread(bdm_res_w8_ZT1)
ZT1<-birdman_dat(ZT1,annotations,"just_annot")
write.table(ZT1,paste0(bdm_dir,"birdman_justcred_wk8ZT1_nonNASHNASH_results_justannot.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)
birdman_plot(ZT1) + labs(title="ZT1 Wk8\nnon-NASH vs. NASH")
ggsave(paste0(bdm_dir,"birdman_justcred_Wk8ZT1_nonNASHNASH_justannot.pdf"),height=6, width=22)

#Wk8 ZT13 204
ZT13<-fread(bdm_res_w8_ZT13)
ZT13<-birdman_dat(ZT13,annotations,"all")
write.table(ZT13,paste0(bdm_dir,"birdman_justcred_wk8ZT13_nonNASHNASH_results.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 
birdman_plot(ZT13) + labs(title="ZT13 Wk8\nnon-NASH vs. NASH")
ggsave(paste0(bdm_dir,"birdman_justcred_Wk8ZT13_nonNASHNASH.pdf"),height=26, width=22)

ZT13<-fread(bdm_res_w8_ZT13)
ZT13<-birdman_dat(ZT13,annotations,"just_annot")
write.table(ZT13,paste0(bdm_dir,"birdman_justcred_wk8ZT13_nonNASHNASH_results_justannot.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 
birdman_plot(ZT13) + labs(title="ZT13 Wk8\nnon-NASH vs. NASH")
ggsave(paste0(bdm_dir,"birdman_justcred_Wk8ZT13_nonNASHNASH_justannot.pdf"),height=8, width=22)

#Wk12 ZT1 466
ZT1<-fread(bdm_res_w12_ZT1)
ZT1<-birdman_dat(ZT1,annotations,"all")
write.table(ZT1,paste0(bdm_dir,"birdman_justcred_wk12ZT1_nonNASHNASH_results.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 
birdman_plot(ZT1) + labs(title="ZT1 Wk12\nnon-NASH vs. NASH")
ggsave(paste0(bdm_dir,"birdman_justcred_Wk12ZT1_nonNASHNASH.pdf"),height=25, width=21)

ZT1<-fread(bdm_res_w12_ZT1)
ZT1<-birdman_dat(ZT1,annotations,"just_annot")
write.table(ZT1,paste0(bdm_dir,"birdman_justcred_wk12ZT1_nonNASHNASH_results_justannot.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)
birdman_plot(ZT1) + labs(title="ZT1 Wk12\nnon-NASH vs. NASH")
ggsave(paste0(bdm_dir,"birdman_justcred_Wk12ZT1_nonNASHNASH_justannot.pdf"),height=12, width=22)

#Wk12 ZT13 634
ZT13<-fread(bdm_res_w12_ZT13)
ZT13<-birdman_dat(ZT13,annotations,"all")
write.table(ZT13,paste0(bdm_dir,"birdman_justcred_wk12ZT13_nonNASHNASH_results.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 
birdman_plot(ZT13) + labs(title="ZT13 Wk12\nnon-NASH vs. NASH")
ggsave(paste0(bdm_dir,"birdman_justcred_Wk12ZT13_nonNASHNASH.pdf"),height=25, width=21)

ZT13<-fread(bdm_res_w12_ZT13)
ZT13<-birdman_dat(ZT13,annotations,"just_annot")
write.table(ZT13,paste0(bdm_dir,"birdman_justcred_wk12ZT13_nonNASHNASH_results_justannot.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)
birdman_plot(ZT13) + labs(title="ZT13 Wk12\nnon-NASH vs. NASH")
ggsave(paste0(bdm_dir,"birdman_justcred_Wk12ZT13_nonNASHNASH_justannot.pdf"),height=18, width=20)

#########################################################################
#Overlap of BIRDMAn results for NASH vs. non-NASH--Fig. 4A

#ZT1 Wk8 to Wk12 (left)

wk8ZT1<-fread(paste0(bdm_dir,"birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv"))
wk12ZT1<-fread(paste0(bdm_dir,"birdman_justcred_wk12ZT1_nonNASHNASH_results.tsv"))

list_venn <- list(Wk8 = wk8ZT1$Name,
                  Wk12 = wk12ZT1$Name)

p<-venn.diagram(list_venn, 
                main="ZT1 Non-NASH vs. NASH",filename = NULL)

grid.draw(p)
pdf(file=paste0(bdm_dir,"birdman_justcredoverlap_nonNASHNASH_ZT1Wk8Wk12results.pdf"))
grid.draw(p)
dev.off()


##ZT13 Wk8 to Wk12 (right)

wk8ZT13<-fread(paste0(bdm_dir,"birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv"))
wk12ZT13<-fread(paste0(bdm_dir,"birdman_justcred_wk12ZT13_nonNASHNASH_results.tsv"))

list_venn <- list(Wk8 = wk8ZT13$Name,
                  Wk12 = wk12ZT13$Name)

p<-venn.diagram(list_venn, 
                main="ZT13 non-NASH vs. NASH",filename = NULL)

grid.draw(p)
pdf(file=paste0(bdm_dir,"birdman_justcredoverlap_nonNASHNASH_ZT13Wk8Wk12results.pdf"))
grid.draw(p)
dev.off()

#########################################################################
#Overlapping plots of BIRDMAn results for NASH vs. non-NASH for ZT1 and ZT13--Fig. 4B

wk8ZT1<-fread(paste0(bdm_dir,"birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv"))
wk12ZT1<-fread(paste0(bdm_dir,"birdman_justcred_wk12ZT1_nonNASHNASH_results.tsv"))
list_venn <- list(Wk8 = wk8ZT1$Name,
                  Wk12 = wk12ZT1$Name)
ItemsList <- venn(list_venn, show.plot = FALSE)
all1<-attributes(ItemsList)$intersections

wk8ZT13<-fread(paste0(bdm_dir,"birdman_justcred_wk8ZT1_nonNASHNASH_results.tsv"))
wk12ZT13<-fread(paste0(bdm_dir,"birdman_justcred_wk12ZT13_nonNASHNASH_results.tsv"))

list_venn <- list(Wk8 = wk8ZT13$Name,
                  Wk12 = wk12ZT13$Name)
ItemsList <- venn(list_venn, show.plot = FALSE)
all2<-attributes(ItemsList)$intersections

#ZT1 (top)
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
ggsave(paste0(bdm_dir,"birdman_justcred_ZT1_wk812comb_shared29_nonNASHNASH.pdf"),height=6 , width=24)
#ggsave(paste0(bdm_dir,"birdman_justcred_ZT1_wk812comb_shared29_ordwk8_nonNASHNASH.pdf"),height=6 , width=24)

#ZT13 (bottom)
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
ggsave(paste0(bdm_dir,"birdman_justcred_ZT13_wk812comb_shared15_nonNASHNASH.pdf"),height=3.5 , width=23)
#ggsave(paste0(bdm_dir,"birdman_justcred_ZT13_wk812comb_shared15_ordwk8_nonNASHNASH.pdf"),height=3.5 , width=22)

#########################################################################
#Normalized abundance plots for metabolites of interest--Fig.4C-D; Fig.S3B-C

mtb<-fread(feattab)%>%
  dplyr::select(c(1,4:5,3,11,21:23,27,33))%>%
  mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition),
         cln_Compound_Name=gsub("\\(predicted.*","",Compound_Name))%>%
  mutate(cln_Compound_Name=gsub("\\(delta.*","",cln_Compound_Name))%>%
  mutate(label_name=ifelse(is.na(cln_Compound_Name),paste("mtb", FeatureID, "m/z",mz, sep=" "),paste("mtb",FeatureID,"m/z",mz, cln_Compound_Name, sep=" ")))%>%
  mutate(ATTRIBUTE_disease_stage=factor(ATTRIBUTE_disease_stage,levels=c("TRF","NASH")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","not_applicable")),
         ATTRIBUTE_condition=factor(ATTRIBUTE_condition,levels=c("NA","FA","FT")))

mtb1264<-paired_plots_ind(mtb,1264,"ZT1")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT1_wk812comb_mtb1264_nonNASHNASH.pdf"),height=4 , width=4)

mtb3597<-paired_plots_ind(mtb,3597,"ZT1")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT1_wk812comb_mtb3597_nonNASHNASH.pdf"),height=4 , width=4)

mtb1992<-paired_plots_ind(mtb,1992,"ZT1")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT1_wk812comb_mtb1992_nonNASHNASH.pdf"),height=4 , width=4)

mtb3748<-paired_plots_ind(mtb,3748,"ZT1")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT1_wk812comb_mtb3748_nonNASHNASH.pdf"),height=4 , width=4)

mtb3552<-paired_plots_ind(mtb,3552,"ZT1")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT1_wk812comb_mtb3552_nonNASHNASH.pdf"),height=4 , width=4)

mtb2245<-paired_plots_ind(mtb,2245,"ZT1")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT1_wk812comb_mtb2245_nonNASHNASH.pdf"),height=4 , width=4)

mtb2568<-paired_plots_ind(mtb,2568,"ZT1")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT1_wk812comb_mtb2568_nonNASHNASH.pdf"),height=4 , width=4)

mtb2691<-paired_plots_ind(mtb,2691,"ZT1")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT1_wk812comb_mtb2691_nonNASHNASH.pdf"),height=4 , width=4)

mtb4446<-paired_plots_ind(mtb,4446,"ZT1")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT1_wk812comb_mtb4446_nonNASHNASH.pdf"),height=4 , width=4)

#ZT13
mtb2245<-paired_plots_ind(mtb,2245,"ZT13")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT13_wk812comb_mtb2245_nonNASHNASH.pdf"),height=4 , width=4)

mtb4168<-paired_plots_ind(mtb,4168,"ZT13")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT13_wk812comb_mtb4168_nonNASHNASH.pdf"),height=4 , width=4)

mtb4200<-paired_plots_ind(mtb,4200,"ZT13")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT13_wk812comb_mtb4200_nonNASHNASH.pdf"),height=4 , width=4)

mtb4166<-paired_plots_ind(mtb,4166,"ZT13")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT13_wk812comb_mtb4166_nonNASHNASH.pdf"),height=4 , width=4)

mtb1522<-paired_plots_ind(mtb,1522,"ZT13")
ggsave(paste0(bdm_dir,"birdman_justcred_ZT13_wk812comb_mtb1522_nonNASHNASH.pdf"),height=4 , width=4)
