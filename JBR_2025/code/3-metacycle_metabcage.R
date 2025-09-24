setwd("~/Notebooks/sfloresr/MASH-TRF/JBR_2025")

library(tidyverse)
library(data.table)
library(MetaCycle)

###########################################################
#inputs
metabdata48h_df<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/first48h_metabcagetable.txt"
dat_filepath<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/data/"
res_filepath<-"~/Notebooks/sfloresr/MASH-TRF/JBR_2025/results/"
###########################################################
#clean files to run metacycle

f48h_ZT<-fread(metabdata48h_df)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         sample_name=paste(condition,zt_time_binned,Animal,sep="_"))

md<-f48h_ZT%>%
  dplyr::select(sample_name,condition,zt_time_binned,Animal)
write.table(md,paste0(dat_filepath,"metabcage_md_metacycle.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

mdNA<-md%>%arrange(zt_time_binned)%>%
  filter(condition=="NA")

f48NA<-f48h_ZT%>%filter(condition=="NA")%>%dplyr::select(-(1:3))%>%
  gather(metabcage,values,-sample_name)%>%
  spread(sample_name,values)%>%
  dplyr::select(metabcage, all_of(mdNA$sample_name))

write.table(f48NA,paste0(dat_filepath,"first48h_metabcagetable_NA_metacycle.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

mdFA<-md%>%arrange(zt_time_binned)%>%
  filter(condition=="FA")
f48FA<-f48h_ZT%>%filter(condition=="FA")%>%dplyr::select(-(1:3))%>%
  gather(metabcage,values,-sample_name)%>%
  spread(sample_name,values)%>%
  dplyr::select(metabcage, all_of(mdFA$sample_name))

write.table(f48FA,paste0(dat_filepath,"first48h_metabcagetable_FA_metacycle.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

mdFT<-md%>%arrange(zt_time_binned)%>%
  filter(condition=="FT")
f48FT<-f48h_ZT%>%filter(condition=="FT")%>%dplyr::select(-(1:3))%>%
  gather(metabcage,values,-sample_name)%>%
  spread(sample_name,values)%>%
  dplyr::select(metabcage, all_of(mdFT$sample_name))

write.table(f48FT,paste0(dat_filepath,"first48h_metabcagetable_FT_metacycle.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

###########################################################
#Metacycle

#on NA
meta2d(infile=paste0(dat_filepath,"first48h_metabcagetable_NA_metacycle.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(res_filepath,"first48h_metabcagetable_NA_metacycle/"),timepoints=rep(seq(1, 21, by=4),times=c(3,3,3,3,3,3)),minper=20,maxper=24)

#on FA
meta2d(infile=paste0(dat_filepath,"first48h_metabcagetable_FA_metacycle.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(res_filepath,"first48h_metabcagetable_FA_metacycle/"),timepoints=rep(seq(1, 21, by=4),times=c(3,3,3,3,3,3)),minper=20,maxper=24)

#on FT
meta2d(infile=paste0(dat_filepath,"first48h_metabcagetable_FT_metacycle.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(res_filepath,"first48h_metabcagetable_FT_metacycle/"),timepoints=rep(seq(1, 21, by=4),times=c(8,8,8,8,8,8)),minper=20,maxper=24)

###########################################################

NA_metacyc<-fread(paste0(res_filepath,"first48h_metabcagetable_NA_metacycle/meta2d_first48h_metabcagetable_NA_metacycle.txt"))%>%
  dplyr::rename(metabcage=CycID)

FA_metacyc<-fread(paste0(res_filepath,"first48h_metabcagetable_FA_metacycle/meta2d_first48h_metabcagetable_FA_metacycle.txt"))%>%
  dplyr::rename(metabcage=CycID)

FT_metacyc<-fread(paste0(res_filepath,"first48h_metabcagetable_FT_metacycle/meta2d_first48h_metabcagetable_FT_metacycle.txt"))%>%
  dplyr::rename(metabcage=CycID)

#RER not cycling in any condition
#Sleep cycling under NA and FT but not FA
#motor activity cycling under NA and FT but not FA 
#mean energy expenditure all cycling
#food consumption (kcal) NA and FA not cycling FT cycling
