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

#################################################################

mtb_human<-fread("NAFLD_human_val/massql_helena/combined_tables/Stephany_matches_ReDU_merging_humans.tsv")%>%
  filter(UBERONBodyPartName=="feces")%>%unique()

sel_mtbs<-c(620,851,933,1025,1264,1425,1522,1846,1974,1992,2245,2444,2448,2467,2568,2658,2691,2762,2793,2931,2934,2938,2944,
            2999,3191,3280,3372,3552,3597,3624,3748,3870,3885,4005,4166,4168,4200,4386,4408,4419,4446)

annot<-fread("mtb/stool/data_files/20231108_gnps_run/DB_result/mz_annot_v2.txt")
#################################################################

#NAFLD study MSV000082374
NAFLD<-mtb_human%>%filter(ATTRIBUTE_DatasetAccession=="MSV000082374")
#unique(NAFLD$Scan) #1264 1425 3552

md_sel<-fread("NAFLD_human_val/caussy_nafld_16s//NAFLD_metadata.tsv")%>%
  filter(sample_name %in% NAFLD$SubjectIdentifierAsRecorded)

NAFLD<-NAFLD%>%
  dplyr::rename(sample_name=SubjectIdentifierAsRecorded)%>%
  left_join(.,md_sel,by="sample_name")

mtb1264<-NAFLD%>%filter(Scan==1264) #m/z 293.2475859 
#6325 in nafld mtb 10E,12Z−octadecadienoic acid
mtb3552<-NAFLD%>%filter(Scan==3552) #m/z 293.2476065
#6325 in nafld mtb 10E,12Z−octadecadienoic acid
mtb1425<-NAFLD%>%filter(Scan==1425) #m/z 407.2398494
#1659 and 1546 in nafld mtb citreviridin 

#quanttable MAFLD mtb
md<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/metadata/clean_fecal_metadata.tsv")
mtb_nafld<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/quantification_table/quantification_table-00000.csv")%>%
  dplyr::select(-c(3:13))%>%
  dplyr::rename(FeatureID=`row ID`)%>%
  gather(filename,peak_area,-FeatureID, -`row m/z`)%>%
  mutate(filename=gsub("\\.mzXML.*", ".mzXML", filename))%>%
  #filter(FeatureID==6325)%>%
  #filter(FeatureID==1659)%>%
  filter(FeatureID==1546)%>%
  left_join(.,md,by="filename")%>%
  filter(ATTRIBUTE_groups=="G1P"|ATTRIBUTE_groups=="G2P" |ATTRIBUTE_groups=="G3P")

#birdman
#bdm<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/ml_analysis/birdman/fecal-ft-matched-filt10.beta_var.tsv")
#6325 is increased in T1 (advance fibrosis) but the range is too wide so not credible
#1659 is increased in T1 (advance fibrosis) but the range is too wide so not credible
#1546 is increased in T1 (advance fibrosis) and the range is too wide so not credible

#boxplot
p <- ggplot(mtb_nafld, aes(x=ATTRIBUTE_groups, y=log10(peak_area+1), fill=ATTRIBUTE_groups)) + 
  geom_violin(trim=FALSE, alpha=0.5)+geom_jitter(shape=16, position=position_jitter(0.2)) +theme_minimal()+
  theme(legend.position = "top")+
  scale_fill_manual(values=c("#0000a7","forestgreen","#c1272d"))

#################################################################

#filter mtb for G1P and G3P to run bdm
md<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/metadata/clean_fecal_metadata.tsv")%>%
  mutate(sample_name=paste(filename," Peak area",sep=""))%>%
  dplyr::select(sample_name,everything())
write.table(md,"/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/metadata/clean_fecal_metadata_forbdm.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

mdg1pg3p<-md%>%filter(ATTRIBUTE_groups=="G1P"|ATTRIBUTE_groups=="G3P")
mdg1pg2p<-md%>%filter(ATTRIBUTE_groups=="G1P"|ATTRIBUTE_groups=="G2P")

mtb_nafld<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/quantification_table/quantification_table-00000.csv")

mtb_nafld_g3p<-mtb_nafld%>%dplyr::select(`row ID`, any_of(mdg1pg3p$sample_name))
write.table(mtb_nafld_g3p,"/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/quantification_table/quantification_table-G1PG3P.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtb_nafld_g2p<-mtb_nafld%>%dplyr::select(`row ID`, any_of(mdg1pg2p$sample_name))
write.table(mtb_nafld_g2p,"/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/quantification_table/quantification_table-G1PG2P.txt",sep = "\t",row.names = FALSE, quote=FALSE)
#################################################################
#bdm results
annot_FA<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/DB_analogresult/e902b4b271a04788b8bc9a1f25b0afd3.tsv")%>%
  filter(grepl("Fatty acids",npclassifier_pathway))

annot_BA<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/DB_analogresult/e902b4b271a04788b8bc9a1f25b0afd3.tsv")%>%
  filter(Organism=="BILELIB19"|Organism=="GNPS-BILE-ACID-MODIFICATIONS")
  
annot<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/metadata/mtb_annotation_key.tsv")%>%
  mutate(Compound_Name = gsub("\\(predicted.*", "", Compound_Name))
mz<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/quantification_table/quantification_table-00000.csv")%>%
  dplyr::select(1:2)%>%dplyr::rename(FeatureID=`row ID`,mz=`row m/z`)

bdm<-fread("NAFLD_human_val/caussy_nafld_mtb/quantification_table-G1PG3P.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G3P]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G3P]_hdi`=gsub("[(]|[)]","",`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G3P]_hdi`))%>%
  separate(`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G3P]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  #cred=ifelse(min>0.5|max< -0.5,"credible","not_credible"))%>%
  #cred=ifelse(min>1|max< -1,"credible","not_credible"))%>%
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

mz<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/quantification_table/quantification_table-00000.csv")%>%
  dplyr::select(1:3)%>%dplyr::rename(FeatureID=`row ID`,mz=`row m/z`, RT=`row retention time`)

annot<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/metadata/mtb_annotation_key.tsv")%>%
  right_join(.,mz,by="FeatureID")

annotations<-fread("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/NAFLD_ML/data_files/FBMN_mzmine3_wBA/DB_analogresult/e902b4b271a04788b8bc9a1f25b0afd3.tsv")%>%
  dplyr::rename(FeatureID=`#Scan#`)%>%
  dplyr::select(FeatureID,Precursor_MZ,RT_Query,Compound_Name,LibraryQualityString,SharedPeaks,LibraryName,
                MQScore,MassDiff,MZErrorPPM,Adduct,IonMode)

fullres<-fread("NAFLD_human_val/caussy_nafld_mtb/quantification_table-G1PG2P.beta_var.tsv")
bdm<-birdman_dat2(fullres,annot,annotations)%>%
  mutate(cat=ifelse(ratio<0,"Non-MAFLD","MAFLD no AF"))%>%
  mutate(cat=factor(cat,levels=c("Non-MAFLD","MAFLD no AF")))
write.table(bdm,"NAFLD_human_val/caussy_nafld_mtb/bdm_G1PG2P_mtbhits_all.txt",sep = "\t",row.names = FALSE, quote=FALSE)


bdm<-fread("NAFLD_human_val/caussy_nafld_mtb/quantification_table-G1PG2P.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_hdi`=gsub("[(]|[)]","",`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_hdi`))%>%
  separate(`C(ATTRIBUTE_groups, Treatment('G1P'))[T.G2P]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  #cred=ifelse(min>0.5|max< -0.5,"credible","not_credible"))%>%
  #cred=ifelse(min>1|max< -1,"credible","not_credible"))%>%
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


