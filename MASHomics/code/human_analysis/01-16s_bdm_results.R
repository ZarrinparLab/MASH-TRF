setwd("~/Notebooks/sfloresr/MASH-TRF/MASHomics/")

library(tidyverse)
library(data.table)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library("qiime2R")

#16S validation on human NAFLD dataset 
##################################################################
#inputs
metadata<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/human_analysis/caussy_nafld_16s/NAFLD_metadata.tsv"
feattab<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/human_analysis/caussy_nafld_16s/04.taxonomy_filtered.asv.counts.qza"
taxonomy<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/human_analysis/caussy_nafld_16s/04.sklearn.gg2.asv.taxonomy.qza"
data_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/human_analysis/caussy_nafld_16s/"
bdm_results_g2p<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/human_analysis/caussy_nafld_16s/birdman/04.taxonomy_filtered.asv.counts.g1pg2p.beta_var.tsv"
bdm_results_g3p<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/human_analysis/caussy_nafld_16s/birdman/04.taxonomy_filtered.asv.counts.g1pg3p.beta_var.tsv"
results_dir<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/human_analysis/caussy_nafld_16s/birdman/"
##################################################################
#cleaning files for BIRDMAn

taxa<-read_qza(taxonomy)$data%>%
  dplyr::rename(FeatureID=Feature.ID)
write.table(taxa,paste0(data_dir,"taxonomy.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

#g1p vs. g2p
md<-fread(metadata)%>%
  filter(groups=="g1p"|groups=="g2p")%>% #NAFLD w/o AF
  mutate(sample_name = sub(".*?\\.", "", sample_name))

set.seed(1234)
data_mod <- md %>% group_by(groups) %>% sample_frac(.3)

md<-md%>%mutate(traintest=ifelse(sample_name %in% data_mod$sample_name, "test","train"))
write.table(md,paste0(data_dir,"NAFLD_metadata_clean_g2p.tsv"),sep = "\t",row.names = FALSE,quote=FALSE)

m16s<-read_qza(feattab)$data%>%
  as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  dplyr::select(FeatureID,all_of(md$sample_name))
write.table(m16s,paste0(data_dir,"04.taxonomy_filtered.asv.counts.g1pg2p.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

#g1p vs. g3p
md<-fread(metadata)%>%
  filter(groups=="g1p"|groups=="g3p")%>% #NAFLD-cirrhosis
  mutate(sample_name = sub(".*?\\.", "", sample_name))

set.seed(1234)
data_mod <- md %>% group_by(groups) %>% sample_frac(.3)

md<-md%>%mutate(traintest=ifelse(sample_name %in% data_mod$sample_name, "test","train"))
write.table(md,paste0(data_dir,"NAFLD_metadata_clean.tsv"),sep = "\t",row.names = FALSE,quote=FALSE)

m16s<-read_qza(feattab)$data%>%
  as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  dplyr::select(FeatureID,all_of(md$sample_name))

write.table(m16s,paste0(data_dir,"04.taxonomy_filtered.asv.counts.g1pg3p.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

##################################################################
#BIRDMAn results g1p vs g2p--Fig. 5B and Table S3
#BDM params
# beta_prior: float = 5.0
# inv_disp_sd: float = 0.5

bdm<-fread(bdm_results_g2p)%>%
  dplyr::rename(ratio=`C(groups, Treatment('g1p'))[T.g2p]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(groups, Treatment('g1p'))[T.g2p]_hdi`=gsub("[(]|[)]","",`C(groups, Treatment('g1p'))[T.g2p]_hdi`))%>%
  separate(`C(groups, Treatment('g1p'))[T.g2p]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>1|max< -1,"credible","not_credible"))%>%
  left_join(.,taxa, by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon, sep=" "))%>%
  dplyr::select(FeatureID, Name,ratio, min, max, cred)%>%
  arrange(ratio)%>%
  mutate(cat=ifelse(ratio<0,"Non-MAFLD","MAFLD no AF"))%>%
  mutate(cat=factor(cat,levels=c("Non-MAFLD","MAFLD no AF")))
write.table(bdm,paste0(results_dir,"birdman_g1pg2p_results.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) #Table S3


bdm<-fread(bdm_results_g2p)%>%
  dplyr::rename(ratio=`C(groups, Treatment('g1p'))[T.g2p]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(groups, Treatment('g1p'))[T.g2p]_hdi`=gsub("[(]|[)]","",`C(groups, Treatment('g1p'))[T.g2p]_hdi`))%>%
  separate(`C(groups, Treatment('g1p'))[T.g2p]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>1|max< -1,"credible","not_credible"))%>%
  left_join(.,taxa, by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon, sep=" "))%>%
  filter(cred=="credible")%>%
  dplyr::select(FeatureID, ratio, min, max, Name)%>%
  arrange(ratio)%>%
  mutate(cat=ifelse(ratio<0,"Non-MAFLD","MAFLD no AF"))%>%
  mutate(cat=factor(cat,levels=c("Non-MAFLD","MAFLD no AF")))
write.table(bdm,paste0(results_dir,"birdman_justcred_g1pg2p_results_cred1.tsv"),sep = "\t",row.names = FALSE, quote=FALSE) 

bdm$Name <- factor(bdm$Name,levels = bdm$Name)

ggplot(bdm, aes(x =Name , y = ratio, ymin = min, ymax = max, fill=cat)) + 
  geom_bar(stat="identity", alpha=0.7)+
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(title="Non-MAFLD vs. MAFLD")+
  coord_flip()
ggsave(paste0(results_dir,"SFR24_0213_g1pg2p_cred1_colored.pdf"),height=4, width=26) #Fig. 5B

##################################################################
#BIRDMAn results g1p vs. g3p--Fig. S4B (top)

bdm<-fread(bdm_results_g3p)%>%
  dplyr::rename(ratio=`C(groups, Treatment('g1p'))[T.g3p]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(groups, Treatment('g1p'))[T.g3p]_hdi`=gsub("[(]|[)]","",`C(groups, Treatment('g1p'))[T.g3p]_hdi`))%>%
  separate(`C(groups, Treatment('g1p'))[T.g3p]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  left_join(.,taxa, by="FeatureID")%>%
  mutate(Name=paste(FeatureID,Taxon, sep=" "))%>%
  filter(cred=="credible")%>%
  dplyr::select(FeatureID, ratio, min, max, Name)%>%
  arrange(ratio)%>%
  mutate(cat=ifelse(ratio<0,"Non-MAFLD","MAFLD"))%>%
  mutate(cat=factor(cat,levels=c("Non-MAFLD","MAFLD")))

write.table(bdm,paste0(results_dir,"birdman_justcred_g1pg3p_results.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

#is there a hit for the lachno or oscillo we saw from our gut microbes here?

blast_oscillo<-fread(paste0(data_dir,"WUMWNPBE114-Alignment_oscillo.txt"))
blast_lachno<-fread(paste0(data_dir,"WUN4T7BF114-Alignment_lachno.txt"))

bdm_oscillo<-bdm%>%
  filter(FeatureID %in% blast_oscillo$`subject acc.ver`)%>%
  mutate(mash_hit="oscillo")

bdm_lachno<-bdm%>%
  filter(FeatureID %in% blast_lachno$`subject acc.ver`)%>%
  mutate(mash_hit="lachno")

combhits<-rbind(bdm_oscillo,bdm_lachno)

ggplot(combhits, aes(x =Name , y = ratio, ymin = min, ymax = max, fill=cat)) + 
  geom_bar(stat="identity", alpha=0.7)+
  geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_fill_manual(values=c("#0000a7","#c1272d"))+
  labs(title="Non-MAFLD vs. MAFLD")+
  coord_flip()
ggsave(paste0(results_dir,"SFR24_0214_g1pg3p_cred0_matchtomouselachnooscillo_colored.pdf"),height=3.5, width=25.5) #Fig. S4B (top)
