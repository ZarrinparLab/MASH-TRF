setwd("~/Notebooks/sfloresr/MASH-TRF/MASHomics/")

library(tidyverse)
library(data.table)
library(ggpubr)
library(VennDiagram)
library("qiime2R")
library(ggalluvial)

##################################################################
#inputs
mash_metadata<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_16s/metadata_cln_addmoreNASHcat.txt"
taxonomy<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/mashstool16s_preprocessed_20211020_ID_13785_gg2/04.sklearn.gg2.asv.taxonomy.qza"
feattab<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/mashstool16s_preprocessed_20211020_ID_13785_gg2/NASH_taxonomy_filtered.gg2.asv.counts-rclr.txt"
data_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/stool_16s/"
results_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/stool_16s/"
##################################################################
#functions

birdman_dat<-function(dat,annot,subdat){
  df<-dat%>%
    dplyr::rename(ratio=`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_mean`,
                FeatureID=Feature)%>%
    mutate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`=gsub("[(]|[)]","",`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`))%>%
    separate(`C(NASH_category, Treatment('Non_NASH'))[T.NASH]_hdi`,c("min","max"), sep=",")%>%
    mutate(min=as.numeric(min),
           max=as.numeric(max),
           cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
    left_join(.,annot, by="FeatureID")%>%
    mutate(Name=paste(FeatureID,Taxon, sep=" "))
  
  if(suddat=="just_credible"){
    df<-df%>%
      filter(cred=="credible")%>%
      dplyr::select(FeatureID, ratio, min, max, Name)%>%
      arrange(ratio)
  }
  else{
    df<-df%>%
      dplyr::select(Name, ratio, min, max, cred)%>%
      arrange(ratio)
  }
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

paired_plots_ind<-function(dat,microid,subdat){
  
  if(subdat=="FA"){
    count_dat<-dat%>%
      filter(FeatureID==microid & condition=="FA" & 
               NASH_category!="not_applicable")
  }
  else{
    count_dat<-dat%>%
      filter(FeatureID==microid & condition=="FT" & 
               NASH_category!="not_applicable")
  }
  
  count_dat<-count_dat%>%
    arrange(collection_timepoint,collection_disease_stage)%>%
    arrange(host_subject_id)%>%as.data.table()
  
  p<-ggplot(count_dat, aes(x =NASH_category , y = counts)) +
    geom_boxplot(aes(fill = NASH_category), alpha = .5) +
    ggtitle(stringr::str_wrap(unique(count_dat$Name), width=30))+
    geom_dotplot(binaxis='y', stackdir='center', dotsize=2)+
    theme_bw()+ theme(legend.position = "top")+
    scale_fill_manual(values=c("#0000a7","#c1272d"))+
    facet_wrap(~ collection_disease_stage)
  return(p)
}

##################################################################
#Combined BIRDMAn results--Table S1

FT<-fread(paste0(results_dir,"birdman_results/NASH_TRF_FT_taxonomy_filtered.gg2.asv.counts.beta_var.tsv"))
FT<-birdman_dat(FT,annot,"all")%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))%>%
  mutate(condition="FT", week=4)

FA<-fread(paste0(results_dir,"birdman_results/NASH_TRF_FA_taxonomy_filtered.gg2.asv.counts.beta_var.tsv"))
FA<-birdman_dat(FA,annot,"all")%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))%>%
  mutate(condition="FA", week=4)

FT2<-fread(paste0(results_dir,"birdman_results/NASH_NASH_FT_taxonomy_filtered.gg2.asv.counts.beta_var.tsv"))
FT2<-birdman_dat(FT2,annot,"all")%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))%>%
  mutate(condition="FT", week=7)

FA2<-fread(paste0(results_dir,"birdman_results/NASH_NASH_FA_taxonomy_filtered.gg2.asv.counts.beta_var.tsv"))
FA2<-birdman_dat(FA2,annot,"all")%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))%>%
  mutate(condition="FA", week=7)

all_res<-rbind(FA,FT,FA2,FT2)
write.table(all_res,paste0(results_dir,"birdman_results/birdman_nonNASHNASH_allresults_gg2.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 

##################################################################
#Differential abundance analysis (BIRDMAn) results NASH vs. non-NASH (FA)--Fig. S2A

annot<-read_qza(taxonomy)$data%>%dplyr::rename(FeatureID=Feature.ID)
write.table(annot,paste0(data_dir,"mashstool16s_preprocessed_20211020_ID_13785_gg2/taxonomy.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 

##FA Wk4 (overall week 8), left
FA<-fread(paste0(results_dir,"birdman_results/NASH_TRF_FA_taxonomy_filtered.gg2.asv.counts.beta_var.tsv"))
FA<-birdman_dat(FA,annot,"just_credible")%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FA,paste0(results_dir,"birdman_results/birdman_justcred_wk8FA_nonNASHNASH_results_gg2.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FA) + labs(title="FA Wk8\nnon-NASH vs. NASH")
ggsave(paste0(results_dir,"birdman_results/birdman_justcred_Wk8FA_nonNASHNASH_gg2.pdf"),height=6, width=14)

##FA Wk7 (overall week 12), right
FA<-fread(paste0(results_dir,"birdman_results/NASH_NASH_FA_taxonomy_filtered.gg2.asv.counts.beta_var.tsv"))
FA<-birdman_dat(FA,annot,"just_credible")%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FA,paste0(results_dir,"birdman_results/birdman_justcred_wk12FA_nonNASHNASH_results_gg2.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FA) + labs(title="FA Wk12\nnon-NASH vs. NASH")
ggsave(paste0(results_dir,"birdman_results/birdman_justcred_Wk12FA_nonNASHNASH_gg2.pdf"),height=7, width=14)

##################################################################
#Differential abundance analysis (BIRDMAn) results NASH vs. non-NASH (FT)--Fig. S2B
##need to run BIRMAn script first

##FT Wk4 (overall week 8), left
FT<-fread(paste0(results_dir,"birdman_results/NASH_TRF_FT_taxonomy_filtered.gg2.asv.counts.beta_var.tsv"))
FT<-birdman_dat(FT,annot,"just_credible")%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FT,paste0(results_dir,"birdman_results/birdman_justcred_wk8FT_nonNASHNASH_results_gg2.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FT) + labs(title="FT Wk8\nnon-NASH vs. NASH")
ggsave(paste0(results_dir,"birdman_results/birdman_justcred_Wk8FT_nonNASHNASH_gg2.pdf"),height=6, width=13)

##FT Wk7 (overall week 12), right
FT<-fread(paste0(results_dir,"birdman_results/NASH_NASH_FT_taxonomy_filtered.gg2.asv.counts.beta_var.tsv"))
FT<-birdman_dat(FT,annot,"just_credible")%>%
  mutate(cat=ifelse(ratio<0,"Non_NASH","NASH"))%>%
  mutate(cat=factor(cat,levels=c("Non_NASH","NASH")))
write.table(FT,paste0(results_dir,"birdman_results/birdman_justcred_wk12FT_nonNASHNASH_results_gg2.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(FT) + labs(title="FT Wk12\nnon-NASH vs. NASH")
ggsave(paste0(results_dir,"birdman_results/birdman_justcred_Wk12FT_nonNASHNASH_gg2.pdf"),height=4, width=13)

##################################################################
#Overlap of BIRDMAn results for NASH vs. non-NASH--Fig. 2A

##FA Wk4 to Wk7, left
wk8FA<-fread(paste0(results_dir,"/birdman_results/birdman_justcred_wk8FA_nonNASHNASH_results_gg2.tsv"))
wk12FA<-fread(paste0(results_dir,"/birdman_results/birdman_justcred_wk12FA_nonNASHNASH_results_gg2.tsv"))

list_venn <- list(Wk8 = wk8FA$Name,
                  Wk12 = wk12FA$Name)

p<-venn.diagram(list_venn, fill = c("#E69F00","#D55E00"),height = 10,lty = 0, 
                main="FA Non-NASH vs. NASH",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file=paste0(results_dir,"birdman_results/birdman_justcredoverlap_nonNASHNASH_FAWk8Wk12results_gg2.pdf"))
grid.draw(p)
dev.off()

##FT Wk8 to Wk12, right
wk8FT<-fread(paste0(results_dir,"birdman_results/birdman_justcred_wk8FT_nonNASHNASH_results_gg2.tsv"))
wk12FT<-fread(paste0(results_dir,"birdman_results/birdman_justcred_wk12FT_nonNASHNASH_results_gg2.tsv"))

list_venn <- list(Wk8 = wk8FT$Name,
                  Wk12 = wk12FT$Name)

p<-venn.diagram(list_venn, fill = c("limegreen","darkgreen"),height = 10,lty = 0, 
                main="FT non-NASH vs. NASH",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file=paste0(results_dir,"birdman_results/birdman_justcredoverlap_nonNASHNASH_FTWk8Wk12results_gg2.pdf"))
grid.draw(p)
dev.off()

##################################################################
#Alluvial plots of BIRDMAn results for NASH vs. non-NASH for FA and FT--Fig. 2B

##FA data for non-nash vs nash going from wk 4 to wk 7, top
FAw12<-fread(paste0(results_dir,"birdman_results/birdman_justcred_wk12FA_nonNASHNASH_results_gg2.tsv"))%>%
  arrange(desc(ratio))%>%rownames_to_column("rank")%>%mutate(week="Wk12",
                                                             Name_rank=paste(rank,Name, sep="."))

FAw8<-fread(paste0(results_dir,"birdman_results/birdman_justcred_wk8FA_nonNASHNASH_results_gg2.tsv"))%>%
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

ggsave(paste0(results_dir,"birdman_results/birdman_justcred_wk8wk12FA_alluvial_gg2.pdf"),height=5, width=8)

##FT data for non-nash vs nash going from wk 4 to wk 7, bottom
FTw12<-fread(paste0(results_dir,"birdman_results/birdman_justcred_wk12FT_nonNASHNASH_results_gg2.tsv"))%>%
  arrange(desc(ratio))%>%rownames_to_column("rank")%>%mutate(week="Wk12",
                                                             Name_rank=paste(rank,Name, sep="."))

FTw8<-fread(paste0(results_dir,"birdman_results/birdman_justcred_wk8FT_nonNASHNASH_results_gg2.tsv"))%>%
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

ggsave(paste0(results_dir,"birdman_results/birdman_justcred_wk8wk12FT_alluvial_gg2.pdf"),height=4, width=8)

##############################################################
#Normalized abundance plots for microbes of interest--Fig.2C-D; Fig.S2C

md<-fread(mash_metadata)
m16s<-fread(feattab)%>%
  dplyr::rename(FeatureID=`#OTU ID`)%>%
  gather(sample_name,counts,-FeatureID)%>%
  mutate(counts=ifelse(is.na(counts),0,counts))%>%
  left_join(.,md,by="sample_name")%>%
  left_join(.,annot,by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon,sep=" "),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH")),
         collection_disease_stage=factor(collection_disease_stage,levels=c("TRF","NASH")))

##FA related hits
micro1<-paired_plots_ind(m16s,"65bc2b52d9d09bde0d63aadfcffa450f","FA")
ggsave(paste0(results_dir,"birdman_results/bdmFA7_65bc2_nonNASHNASH_gg2_unprd.pdf"),height=4 , width=4)

micro2<-paired_plots_ind(m16s,"7f782768e8d3067af0a923af092239d3","FA")
ggsave(paste0(results_dir,"birdman_results/bdmFA7_7f782_nonNASHNASH_gg2_unprd.pdf"),height=4 , width=4)

micro3<-paired_plots_ind(m16s,"726fe9527294f50f4d2654b3a25730a9","FA")
ggsave(paste0(results_dir,"birdman_results/bdmFA7_726fe_nonNASHNASH_gg2_unprd.pdf"),height=4 , width=4)

micro4<-paired_plots_ind(m16s,"eb4104bfe7ff965e7288b9e381317aa8","FA")
ggsave(paste0(results_dir,"birdman_results/bdmFA7_eb410_nonNASHNASH_gg2_unprd.pdf"),height=4 , width=4)

micro5<-paired_plots_ind(m16s,"d1c1529acbbd977e6219701a664424f0","FA")
ggsave(paste0(results_dir,"birdman_results/bdmFA7_d1c15_nonNASHNASH_gg2_unprd.pdf"),height=4 , width=4)

micro6<-paired_plots_ind(m16s,"e613626b2ce8bfbf28995d9b11edaeed","FA")
ggsave(paste0(results_dir,"birdman_results/bdmFA7_e6136_nonNASHNASH_gg2_unprd.pdf"),height=4 , width=4)

micro7<-paired_plots_ind(m16s,"d90a306bb67d75602d9e789894d2c81e","FA")
ggsave(paste0(results_dir,"birdman_results/bdmFA7_d90a3_nonNASHNASH_gg2_unprd.pdf"),height=4 , width=4)

##FT related hits
micro1<-paired_plots_ind(m16s,"d1c1529acbbd977e6219701a664424f0","FT")
ggsave(paste0(results_dir,"birdman_results/bdmFT4_d1c15_nonNASHNASH_gg2_unprd.pdf"),height=4 , width=4)

micro2<-paired_plots_ind(m16s,"a6c1d5c393eaa8b78baddc67c010b5fe","FT")
ggsave(paste0(results_dir,"birdman_results/bdmFT4_a6c1d_nonNASHNASH_gg2_unprd.pdf"),height=4 , width=4)

micro3<-paired_plots_ind(m16s,"be65daef90f69cfc35f5bcd8e387d377","FT")
ggsave(paste0(results_dir,"birdman_results/bdmFT4_be65d_nonNASHNASH_gg2_unprd.pdf"),height=4 , width=4)

micro4<-paired_plots_ind(m16s,"923a477c342aa445c4eecca0f8415ed5","FT")
ggsave(paste0(results_dir,"birdman_results/bdmFT4_923a4_nonNASHNASH_gg2_unprd.pdf"),height=4 , width=4)
