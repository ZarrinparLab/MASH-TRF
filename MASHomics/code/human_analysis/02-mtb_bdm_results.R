setwd("~/Notebooks/sfloresr/MASH-TRF/MASHomics/")

library(tidyverse)
library(data.table)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library("qiime2R")

#mtb validation on human NAFLD dataset--G2P
#################################################################
#inputs
metadata<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/human_analysis/caussy_nafld_mtb/clean_fecal_metadata.tsv"
feattab<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/human_analysis/caussy_nafld_mtb/FBMN_mzmine3_wBA/quantification_table/quantification_table-00000.csv"
annotation<-"/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/metadata/mtb_annotation_key.tsv"
annotation2<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/human_analysis/caussy_nafld_mtb/FBMN_mzmine3_wBA/DB_analogresult/e902b4b271a04788b8bc9a1f25b0afd3.tsv"
data_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/human_analysis/caussy_nafld_mtb/"
bdm_results_g2p<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/human_analysis/caussy_nafld_mtb/quantification_table-G1PG2P.beta_var.tsv"
bdm_results_g3p<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/human_analysis/caussy_nafld_mtb/quantification_table-G1PG3P.beta_var.tsv"
results_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/human_analysis/caussy_nafld_mtb/"
#################################################################
#functions

birdman_dat2<-function(dat,annot,annot2){
  df<-dat%>%
    dplyr::rename(ratio=`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_mean`,
                  FeatureID=Feature)%>%
    mutate(`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_hdi`=gsub("[(]|[)]","",`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_hdi`))%>%
    separate(`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_hdi`,c("min","max"), sep=",")%>%
    mutate(min=as.numeric(min),
           max=as.numeric(max),
           cred=ifelse(min>0|max< 0,"credible","not_credible"))%>%
    left_join(annot,by="FeatureID")%>%
    mutate(label_name=ifelse(is.na(Compound_Name),paste("mtb ",FeatureID," m/z ",signif(mz,5), sept=""),
                             paste("mtb ",FeatureID," m/z ",signif(mz,5), Compound_Name,sept="")))%>%
    left_join(.,annot2,by="FeatureID")%>%
    dplyr::select(FeatureID,label_name,mz,RT,Compound_Name.x,LibraryQualityString,SharedPeaks,LibraryName,
                  MQScore,MassDiff,MZErrorPPM,Adduct,IonMode,ratio, min, max,cred)%>%
    arrange(ratio)
  return(df)
}
#################################################################
#cleaning files
#filter mtb for G1P vs. G2P and G1P vs. G3P to run BIRDMAn

md<-fread(metadata)%>%
  mutate(sample_name=paste(filename," Peak area",sep=""))%>%
  dplyr::select(sample_name,everything())
write.table(md,paste0(data_dir,"metadata/clean_fecal_metadata_forbdm.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

mdg1pg3p<-md%>%filter(ATTRIBUTE_groups=="G1P"|ATTRIBUTE_groups=="G3P")
mdg1pg2p<-md%>%filter(ATTRIBUTE_groups=="G1P"|ATTRIBUTE_groups=="G2P")

mtb_nafld<-fread(feattab)

mtb_nafld_g3p<-mtb_nafld%>%dplyr::select(`row ID`, any_of(mdg1pg3p$sample_name))
write.table(mtb_nafld_g3p,paste0(data_dir,"FBMN_mzmine3_wBA/quantification_table/quantification_table-G1PG3P.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

mtb_nafld_g2p<-mtb_nafld%>%dplyr::select(`row ID`, any_of(mdg1pg2p$sample_name))
write.table(mtb_nafld_g2p,paste0(data_dir,"FBMN_mzmine3_wBA/quantification_table/quantification_table-G1PG2P.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
#################################################################
#BIRDMAn results--Fig. 5C and Table S4

annot<-fread(annotation)%>%
  mutate(Compound_Name = gsub("\\(predicted.*", "", Compound_Name))

mz<-fread(feattab)%>%
  dplyr::select(1:2)%>%dplyr::rename(FeatureID=`row ID`,mz=`row m/z`)

annotations<-fread(annotations2)%>%
  dplyr::rename(FeatureID=`#Scan#`)%>%
  dplyr::select(FeatureID,Precursor_MZ,RT_Query,Compound_Name,LibraryQualityString,SharedPeaks,LibraryName,
                MQScore,MassDiff,MZErrorPPM,Adduct,IonMode)

fullres<-fread(bdm_results_g2p)
bdm<-birdman_dat2(fullres,annot,annotations)%>%
  mutate(cat=ifelse(ratio<0,"Non-MAFLD","MAFLD no AF"))%>%
  mutate(cat=factor(cat,levels=c("Non-MAFLD","MAFLD no AF")))
write.table(bdm,paste0(results_dir,"bdm_G1PG2P_mtbhits_all.txt"),sep = "\t",row.names = FALSE, quote=FALSE) #Table S4

bdm<-fread(bdm_results_g3p)%>%
  dplyr::rename(ratio=`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G3P]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G3P]_hdi`=gsub("[(]|[)]","",`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G3P]_hdi`))%>%
  separate(`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G3P]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  left_join(.,mz, by="FeatureID")%>%
  left_join(annot,by="FeatureID")%>%
  mutate(label_name=ifelse(is.na(Compound_Name),paste("mtb ",FeatureID," m/z ",signif(mz,5), sept=""),
                           paste("mtb ",FeatureID," m/z ",signif(mz,5), Compound_Name,sept="")))%>%
  filter(cred=="credible")%>%
  filter(!is.na(Compound_Name))%>%
  dplyr::select(FeatureID, ratio, min, max, label_name)%>%
  arrange(ratio)%>%
  mutate(cat=ifelse(ratio<0,"Non-MAFLD","MAFLD"))%>%
  mutate(cat=factor(cat,levels=c("Non-MAFLD","MAFLD")))%>%
  filter(FeatureID %in% c(1546,4167,3962))%>%
  mutate(FeatureID=as.character(FeatureID))

write.table(bdm,"NAFLD_human_val/caussy_nafld_mtb/bdm_G1PG3P_cred_mtbhits.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(bdm,"NAFLD_human_val/caussy_nafld_mtb/bdm_G1PG3P_cred_mtbhits_justannot.txt",sep = "\t",row.names = FALSE, quote=FALSE)

bdm<-fread(bdm_results_g2p)%>%
  dplyr::rename(ratio=`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_hdi`=gsub("[(]|[)]","",`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_hdi`))%>%
  separate(`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  left_join(.,mz, by="FeatureID")%>%
  left_join(annot,by="FeatureID")%>%
  mutate(label_name=ifelse(is.na(Compound_Name),paste("mtb ",FeatureID," m/z ",signif(mz,5), sept=""),
                           paste("mtb ",FeatureID," m/z ",signif(mz,5), Compound_Name,sept="")))%>%
  filter(cred=="credible")%>%
  filter(!is.na(Compound_Name))%>%
  dplyr::select(FeatureID, ratio, min, max, label_name)%>%
  arrange(ratio)%>%
  #filter(FeatureID %in% c(1546,4167,3962))%>%
  mutate(cat=ifelse(ratio<0,"Non-MAFLD","MAFLD"))%>%
  mutate(cat=factor(cat,levels=c("Non-MAFLD","MAFLD")))

write.table(bdm,"NAFLD_human_val/caussy_nafld_mtb/bdm_G1PG2P_cred_mtbhits.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(bdm,"NAFLD_human_val/caussy_nafld_mtb/bdm_G1PG2P_cred_mtbhits_justannot.txt",sep = "\t",row.names = FALSE, quote=FALSE)

bdm$label_name <- factor(bdm$label_name,levels =bdm$label_name)
bdm$FeatureID <- factor(bdm$FeatureID,levels =bdm$FeatureID)

# bdm_acids<-bdm%>%filter(grepl("acid",label_name))%>%
#   arrange(ratio)%>%
#   mutate(cat=factor(cat,levels=c("Non-MAFLD","MAFLD")))

#bdm_acids$label_name <- factor(bdm_acids$label_name,levels =bdm_acids$label_name)

p<-ggplot(bdm, aes(x =FeatureID , y = ratio, ymin = min, ymax = max, fill=cat)) + 
  geom_bar(stat="identity", alpha=0.7)+
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  scale_y_continuous(breaks = seq(-15,15, by = 2)) +
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(title="Non-MAFLD vs. MAFLD no AF")+
  coord_flip()
# ggsave("NAFLD_human_val/caussy_nafld_mtb/SFR24_0612_g1pg3p_wsuspectlib.pdf", plot=p,height=8, width=24)
# ggsave("NAFLD_human_val/caussy_nafld_mtb/SFR24_0612_g1pg3p_justbdm_sighits_wmousematch.pdf", plot=p,height=2.5, width=24)
# ggsave("NAFLD_human_val/caussy_nafld_mtb/SFR24_0612_g1pg3p_justbdm_sighits_wmousematch_notannot.pdf", plot=p,height=1.75, width=4)
ggsave("NAFLD_human_val/caussy_nafld_mtb/SFR26_0330_g1pg2p_wsuspectlib.pdf", plot=p,height=8, width=5)

p<-ggplot(bdm, aes(x =label_name , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("Non-MAFLD vs. MAFLD")+
  coord_flip()
ggsave("NAFLD_human_val/caussy_nafld_mtb/SFR24_0513_g1pg2p_BAFA_wsuspectlib.pdf", plot=p,height=6, width=14)
ggsave("NAFLD_human_val/caussy_nafld_mtb/SFR26_0330_g1pg2p_wsuspectlib_annot.pdf", plot=p,height=8, width=14)
ggsave("NAFLD_human_val/caussy_nafld_mtb/SFR26_0330_g1pg2p_wsuspectlib_annot_0.5.pdf", plot=p,height=3.5, width=12.5)

bdm_sel<-bdm%>%filter(FeatureID %in% c(6325,1659,1546))%>%
  arrange(ratio)%>%
  mutate(cat=factor(cat,levels=c("Non-MAFLD","MAFLD")))

bdm_sel$label_name <- factor(bdm_sel$label_name,levels = bdm_sel$label_name)

p<-ggplot(bdm_sel, aes(x =label_name , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("Non-MAFLD vs. MAFLD")+
  coord_flip()
ggsave("NAFLD_human_val/caussy_nafld_mtb/SFR24_0513_g1pg2p_frommMASHstudy.pdf", plot=p,height=2, width=5)

bdm_FA<-bdm%>%filter(FeatureID %in% annot_FA$`#Scan#`)%>%
  mutate(annot="fatty acids")
bdm_BA<-bdm%>%filter(FeatureID %in% annot_BA$`#Scan#`)%>%
  mutate(annot="bile acids")

combFABA<-rbind(bdm_FA,bdm_BA)%>%
  arrange(ratio)%>%
  mutate(cat=factor(cat,levels=c("Non-MAFLD","MAFLD")))

combFABA$label_name <- factor(combFABA$label_name,levels = combFABA$label_name)

p<-ggplot(combFABA, aes(x =label_name , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("Non-MAFLD vs. MAFLD")+
  coord_flip()
ggsave("NAFLD_human_val/caussy_nafld_mtb/SFR24_0513_g1pg2p_BAFA.pdf", plot=p,height=3, width=12)

#################################################################
#results from Helena, matches from human to mouse data

annot_masst<-fread("NAFLD_human_val/caussy_nafld_mtb/masst_human_to_mouse/Annotated_compounds_cos07_tolerance002.txt")%>%
  filter(features!="")%>%
  separate(features,c("mouse_dat","mouse_feat","vs","human_dat","human_feat"),"_")%>%
  dplyr::select(-vs)

notannot_masst<-fread("NAFLD_human_val/caussy_nafld_mtb/masst_human_to_mouse/Not_Annotated_compounds_cos07_tolerance002.txt")%>%
  separate(features,c("mouse_dat","mouse_feat","vs","human_dat","human_feat"),"_")%>%
  dplyr::select(-vs)

masst_results<-rbind(annot_masst,notannot_masst)
write.table(masst_results,"NAFLD_human_val/caussy_nafld_mtb/masst_human_to_mouse/comb_masst_results.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#are there matches to our 41 sel mtbs?

mtb_sel41<-fread("multiomics_16smtb/jrpca_mtbbdm_sel41_feattab.txt")%>%
  filter(FeatureID %in% annot_masst$mouse_feat) #1425, 2691
  #filter(FeatureID %in% notannot_masst$mouse_feat) #none

#1425 goes with 1546
#2691 goes with 4167 and 3962


