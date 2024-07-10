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
library(cowplot)
#library(vegan)
#mtb analysis based on condition
#########################################################################
birdman_dat<-function(dat,annot,compar){
  if(compar=="FAFT"){
    df<-dat%>%
      dplyr::rename(ratio=`C(ATTRIBUTE_condition, Treatment('FA'))[T.FT]_mean`,
                    FeatureID=Feature)%>%
      mutate(`C(ATTRIBUTE_condition, Treatment('FA'))[T.FT]_hdi`=gsub("[(]|[)]","",`C(ATTRIBUTE_condition, Treatment('FA'))[T.FT]_hdi`))%>%
      separate(`C(ATTRIBUTE_condition, Treatment('FA'))[T.FT]_hdi`,c("min","max"), sep=",")
    
  }
  else if(compar=="NAFA"){
    df<-dat%>%
      dplyr::rename(ratio=`C(ATTRIBUTE_condition, Treatment('NA'))[T.FA]_mean`,
                    FeatureID=Feature)%>%
      mutate(`C(ATTRIBUTE_condition, Treatment('NA'))[T.FA]_hdi`=gsub("[(]|[)]","",`C(ATTRIBUTE_condition, Treatment('NA'))[T.FA]_hdi`))%>%
      separate(`C(ATTRIBUTE_condition, Treatment('NA'))[T.FA]_hdi`,c("min","max"), sep=",")  
  }
  else{
    df<-dat%>%
      dplyr::rename(ratio=`C(ATTRIBUTE_condition, Treatment('NA'))[T.FT]_mean`,
                    FeatureID=Feature)%>%
      mutate(`C(ATTRIBUTE_condition, Treatment('NA'))[T.FT]_hdi`=gsub("[(]|[)]","",`C(ATTRIBUTE_condition, Treatment('NA'))[T.FT]_hdi`))%>%
      separate(`C(ATTRIBUTE_condition, Treatment('NA'))[T.FT]_hdi`,c("min","max"), sep=",") 
  }
  df<-df%>%
    mutate(min=as.numeric(min),
           max=as.numeric(max),
           #cred=ifelse(min>0.5|max< -0.5,"credible","not_credible"))%>%
           cred=ifelse(min>0|max< 0,"credible","not_credible"))%>%
    left_join(.,annot, by="FeatureID")%>%
    mutate(Name=paste(FeatureID,Compound_Name,mtb_group,sep=" "))%>%
    mutate(Name=ifelse(is.na(Compound_Name),paste(FeatureID, "m/z",mz,mtb_group,sep=" "),Name))%>%
    filter(cred=="credible") %>%
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
#########################################################################
mz_info<-fread('mtb/stool/data_files/20231108_gnps_run/DB_result/mz_annot_v2.txt')

annotations<-fread("mtb/stool/data_files/20231108_gnps_run/DB_result/DB_results_clean.txt")%>%
  select(FeatureID,Compound_Name)%>%
  mutate(FeatureID=as.integer(FeatureID),
         Compound_Name=gsub("\\(predic.*","",Compound_Name),
         Compound_Name=gsub('"',"",Compound_Name))%>%
  right_join(.,mz_info,by="FeatureID")

#Wk12 FAFT
ZT13FAFT<-fread("mtb/stool/birdman/just_conditions/quantification_table-00000-clean-mash-Wk12-ZT13FAFT_log.beta_var.tsv")
ZT13FAFT<-birdman_dat(ZT13FAFT,annotations,"FAFT")%>%
  mutate(cat=ifelse(ratio<0,"FA","FT"))%>%
  mutate(cat=factor(cat,levels=c("FA","FT")))
write.table(ZT13FAFT,"mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT13_FAFT_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT13FAFT) + labs(title="ZT13 Wk12 FA vs. FT")
ggsave("mtb/stool/birdman/just_conditions/birdman_justcred_Wk12ZT13_FAFT.pdf",height=6, width=13)


ZT1FAFT<-fread("mtb/stool/birdman/just_conditions/quantification_table-00000-clean-mash-Wk12-ZT1FAFT_log.beta_var.tsv")
ZT1FAFT<-birdman_dat(ZT1FAFT,annotations,"FAFT")%>%
  mutate(cat=ifelse(ratio<0,"FA","FT"))%>%
  mutate(cat=factor(cat,levels=c("FA","FT")))
write.table(ZT1FAFT,"mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT1_FAFT_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT1FAFT) + labs(title="ZT1 Wk12 FA vs. FT")
ggsave("mtb/stool/birdman/just_conditions/birdman_justcred_Wk12ZT1_FAFT.pdf",height=13, width=13)

#Wk12 NAFA
ZT13NAFA<-fread("mtb/stool/birdman/just_conditions/quantification_table-00000-clean-mash-Wk12-ZT13NAFA_log.beta_var.tsv")
ZT13NAFA<-birdman_dat(ZT13NAFA,annotations,"NAFA")%>%
  mutate(cat=ifelse(ratio<0,"NA","FA"))%>%
  mutate(cat=factor(cat,levels=c("NA","FA")))%>%
  filter(!grepl("m/z",Name))%>%
  filter(ratio>1 | ratio< -1)
write.table(ZT13NAFA,"mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT13_NAFA_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
write.table(ZT13NAFA,"mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT13_NAFA_results_annot.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT13NAFA) + labs(title="ZT13 Wk12 NA vs. FA")
ggsave("mtb/stool/birdman/just_conditions/birdman_justcred_Wk12ZT13_NAFA_annot_ratio1.pdf",height=10, width=16)


ZT1NAFA<-fread("mtb/stool/birdman/just_conditions/quantification_table-00000-clean-mash-Wk12-ZT1NAFA_log.beta_var.tsv")
ZT1NAFA<-birdman_dat(ZT1NAFA,annotations,"NAFA")%>%
  mutate(cat=ifelse(ratio<0,"NA","FA"))%>%
  mutate(cat=factor(cat,levels=c("NA","FA")))%>%
  filter(!grepl("m/z",Name))%>%
  filter(ratio>1 | ratio< -1)
write.table(ZT1NAFA,"mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT1_NAFA_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
write.table(ZT1NAFA,"mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT1_NAFA_results_annot.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT1NAFA) + labs(title="ZT1 Wk12 NA vs. FA")
ggsave("mtb/stool/birdman/just_conditions/birdman_justcred_Wk12ZT1_NAFA_ratio1.pdf",height=10, width=16)

#Wk12 NAFA
ZT13NAFT<-fread("mtb/stool/birdman/just_conditions/quantification_table-00000-clean-mash-Wk12-ZT13NAFT_log.beta_var.tsv")
ZT13NAFT<-birdman_dat(ZT13NAFT,annotations,"NAFT")%>%
  mutate(cat=ifelse(ratio<0,"NA","FT"))%>%
  mutate(cat=factor(cat,levels=c("NA","FT")))%>%
  filter(!grepl("m/z",Name))%>%
  filter(ratio>1 | ratio< -1)
write.table(ZT13NAFT,"mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT13_NAFT_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
write.table(ZT13NAFT,"mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT13_NAFT_results_annot.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT13NAFT) + labs(title="ZT13 Wk12 NA vs. FT")
ggsave("mtb/stool/birdman/just_conditions/birdman_justcred_Wk12ZT13_NAFT_annot_ratio1.pdf",height=10, width=16)


ZT1NAFT<-fread("mtb/stool/birdman/just_conditions/quantification_table-00000-clean-mash-Wk12-ZT1NAFT_log.beta_var.tsv")
ZT1NAFT<-birdman_dat(ZT1NAFT,annotations,"NAFT")%>%
  mutate(cat=ifelse(ratio<0,"NA","FT"))%>%
  mutate(cat=factor(cat,levels=c("NA","FT")))%>%
  filter(!grepl("m/z",Name))%>%
  filter(ratio>1 | ratio< -1)
write.table(ZT1NAFT,"mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT1_NAFT_results.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 
write.table(ZT1NAFT,"mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT1_NAFT_results_annot.tsv",sep = "\t",row.names = FALSE, quote=FALSE) 

birdman_plot(ZT1NAFT) + labs(title="ZT1 Wk12 NA vs. FT")
ggsave("mtb/stool/birdman/just_conditions/birdman_justcred_Wk12ZT1_NAFT_ratio1.pdf",height=10, width=16)
#########################################################################

#overlap

ZT13_FAFT<-fread('mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT13_FAFT_results.tsv')
ZT13_NAFA<-fread('mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT13_NAFA_results.tsv')
ZT13_NAFT<-fread('mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT13_NAFT_results.tsv')

list_venn <- list(NAFA = ZT13_NAFA$FeatureID,
                  NAFT = ZT13_NAFT$FeatureID,
                  FAFT = ZT13_FAFT$FeatureID)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections
ZT13_FTs<-all$`NAFT:FAFT`
annot_ZT13FT<-annotations%>%filter(FeatureID %in%ZT13_FTs)

p<-venn.diagram(list_venn, height = 10,lty = 1,
                main="ZT13 conditions",
                width = 10,lwd =1,filename = NULL)

grid.draw(p)
pdf(file="mtb/stool/birdman/just_conditions/bdm_venn_ZT13conds.pdf")
grid.draw(p)
dev.off()

ZT1_FAFT<-fread('mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT1_FAFT_results.tsv')
ZT1_NAFA<-fread('mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT1_NAFA_results.tsv')
ZT1_NAFT<-fread('mtb/stool/birdman/just_conditions/birdman_justcred_wk12ZT1_NAFT_results.tsv')

list_venn <- list(NAFA = ZT1_NAFA$FeatureID,
                  NAFT = ZT1_NAFT$FeatureID,
                  FAFT = ZT1_FAFT$FeatureID)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections
ZT1_FTs<-all$`NAFT:FAFT`
annot_ZT1FT<-annotations%>%filter(FeatureID %in%ZT1_FTs)

p<-venn.diagram(list_venn, height = 10,lty = 1,
                main="ZT1 conditions",
                width = 10,lwd =1,filename = NULL)

grid.draw(p)
pdf(file="mtb/stool/birdman/just_conditions/bdm_venn_ZT1conds.pdf")
grid.draw(p)
dev.off()

list_venn <- list(ZT1=ZT1_FAFT$FeatureID,
                  ZT13=ZT13_FAFT$FeatureID)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

p<-venn.diagram(list_venn, height = 10,lty = 1,
                main="FAFT",
                width = 10,lwd =1,filename = NULL)

grid.draw(p)
pdf(file="mtb/stool/birdman/just_conditions/bdm_venn_ZT_FAFT.pdf")
grid.draw(p)
dev.off()

list_venn <- list(ZT1=ZT1_NAFA$FeatureID,
                  ZT13=ZT13_NAFA$FeatureID)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

p<-venn.diagram(list_venn, height = 10,lty = 1,
                main="NAFA",
                width = 10,lwd =1,filename = NULL)

grid.draw(p)
pdf(file="mtb/stool/birdman/just_conditions/bdm_venn_ZT_NAFA.pdf")
grid.draw(p)
dev.off()

list_venn <- list(ZT1=ZT1_NAFT$FeatureID,
                  ZT13=ZT13_NAFT$FeatureID)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

p<-venn.diagram(list_venn, height = 10,lty = 1,
                main="NAFT",
                width = 10,lwd =1,filename = NULL)

grid.draw(p)
pdf(file="mtb/stool/birdman/just_conditions/bdm_venn_ZT_NAFT.pdf")
grid.draw(p)
dev.off()

#########################################################################
md<-fread("mtb/stool/NASH_stool_metabolomics_metadata_new.txt")%>%
  mutate(ATTRIBUTE_condition=ifelse(is.na(ATTRIBUTE_condition),"NA",ATTRIBUTE_condition))%>%
  filter(ATTRIBUTE_disease_stage=="NASH")%>%
  mutate(ATTRIBUTE_condition=factor(ATTRIBUTE_condition,levels=c("NA","FA","FT")))

mtb<-fread("mtb/stool/data_files/20231108_gnps_run/quantification_table/quantification_table-00000-clean-mash-rmblnk.txt")%>%
  dplyr::select(`row ID`, any_of(md$sample_name))%>%
  dplyr::rename(FeatureID=`row ID`)%>%
  gather(sample_name,peak_abun,-FeatureID)%>%
  filter(FeatureID==1184|FeatureID==3471|FeatureID==2697)%>%
  left_join(.,md,by="sample_name")%>%
  left_join(.,annotations,by="FeatureID")

p<-ggplot(mtb, aes(x=ATTRIBUTE_condition, y=log10(peak_abun+1), fill=ATTRIBUTE_timepoint)) +
  geom_boxplot()+ facet_wrap(~FeatureID)

mtb_sub<-mtb%>%
  filter(ATTRIBUTE_timepoint=="ZT1")
p<-ggplot(mtb_sub, aes(x=ATTRIBUTE_condition, y=log10(peak_abun+1), fill=NASH_category)) +
  geom_boxplot()+ facet_wrap(~FeatureID)

