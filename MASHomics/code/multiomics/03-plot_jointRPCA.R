setwd("~/Notebooks/sfloresr/MASH-TRF/MASHomics/")

library(tidyverse)
library(data.table)
library(pheatmap)
library(viridis)

####################################################################
#inputs
annot_16s<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/multiomics/jrpca_16s_annotation_key_gg2.txt"
annot_mtb<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/multiomics/jrpca_mtb_annotation_key_wcanopus.txt"
sel16s_NASHcat_p<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/multiomics/jrpca_sel16s_NASHcat.txt"
sel41_ZTcat_p<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/data/multiomics/jrpca_sel41_ZTcat.txt"
jrpca_res_wk12<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/multiomics/result_bdm_sel41mtb_sel10micro_gg2/"
jrpca_res_wk8<-"~/Notebooks/sfloresr/MASH-TRF/MASHomics/results/multiomics/result_bdm_sel41mtb_sel10micro_wk8_gg2/"
####################################################################
#plot the joint-RPCA correlation (wk12)--Fig. 5A

sel16s_NASHcat<-fread(sel16s_NASHcat_p)
sel41_ZTcat<-fread(sel41_ZTcat_p)

m16s_key<-fread(annot_16s)%>%
  right_join(.,sel16s_NASHcat,by="FeatureID")%>%
  mutate(name=paste(FeatureID,Taxon))%>%
  dplyr::select(FeatureID,name)

mtb_key<-fread(annot_mtb)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  right_join(.,sel41_ZTcat,by="FeatureID")%>%
  mutate(name=ifelse(is.na(Compound_Name),
                     paste("mtb",FeatureID,"m/z",signif(mz,4), sep=" "),
                     paste("mtb",FeatureID,"m/z",signif(mz,4), Compound_Name,sep=" ")),
         FeatureID=as.character(FeatureID))%>%
  dplyr::select(FeatureID,name)

key<-rbind(m16s_key,mtb_key,fill=TRUE)

corr_n<-fread(paste0(jrpca_res_wk12,"correlation_table/Correlation.tsv"),header=TRUE)%>%
  dplyr::rename(FeatureID=featureid)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(53,2:11))%>%
  column_to_rownames("name")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(53,12:52))%>%
  column_to_rownames("name")

write.table(corr_n, paste0(jrpca_res_wk12,"correlation_table/Correlation_clean.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

m16s_key_ext<-fread(annot_16s)%>%
  right_join(.,sel16s_NASHcat,by="FeatureID")%>%
  mutate(name=paste(FeatureID,Taxon),
         NASH_bdmcat=factor(NASH_bdmcat,levels=c("Non-MASH","MASH","Non-MASH:MASH")),
         condition=factor(condition,levels=c("FA","FT","FA:FT")))%>%
  dplyr::select(name,NASH_bdmcat,condition)%>%
  column_to_rownames("name")
mtb_key_ext<-fread(annot_mtb)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  right_join(.,sel41_ZTcat,by="FeatureID")%>%
  mutate(name=ifelse(is.na(Compound_Name),
                     paste("mtb",FeatureID,"m/z",signif(mz,4), sep=" "),
                     paste("mtb",FeatureID,"m/z",signif(mz,4), Compound_Name,sep=" ")),
         FeatureID=as.character(FeatureID),
         mtb_group=factor(mtb_group,levels=c("bile acids","fatty acids & other lipids","unknown")),
         ZT_cat=factor(ZT_cat,levels=c("ZT1","ZT13","ZT1:ZT13")))%>%
  dplyr::select(name,mtb_group,ZT_cat)%>%
  column_to_rownames("name")

ann_colors <- list(
  NASH_bdmcat = c("Non-MASH" = "#0000a7", "MASH" = "#c1272d", "Non-MASH:MASH" = "#92278F"),
  condition = c("FA" = "#D55E00", "FT" = "#009E73", "FA:FT" = "grey50"),
  mtb_group = c("bile acids" = "#006838", "fatty acids & other lipids" = "#F15A29", "unknown" = "#808285"),
  ZT_cat = c("ZT1" = "grey90", "ZT13" = "grey10", "ZT1:ZT13" = "grey50")
)

p<-pheatmap(corr_n, annotation_row = m16s_key_ext,annotation_col = mtb_key_ext,
            annotation_colors = ann_colors,color= inferno(10), main="Wk12")
ggsave(paste0(jrpca_res_wk12,"correlation_table/SFR24_0212_heatmap_wbdmannot.pdf"), p,width = 23, height = 26)

ordw12_micro<-rownames(corr_n[p$tree_row[["order"]],]) #microbe order
cordw12_mtb<-colnames(corr_n[,p$tree_col[["order"]]]) #mtb order
#######################################################################################
#plot the joint-RPCA correlation (wk8)--Fig. S4A

corr_n<-fread(paste0(jrpca_res_wk8,"correlation_table/Correlation.tsv"),header=TRUE)%>%
  dplyr::rename(FeatureID=featureid)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(53,2:11))%>%
  column_to_rownames("name")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(53,12:52))%>%
  column_to_rownames("name")

write.table(corr_n, paste0(jrpca_res_wk8,"correlation_table/Correlation_clean.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

#keep the order of wk12 heatmap (turn off clustering)
corr_n<-corr_n%>%dplyr::select(all_of(cordw12_mtb))%>%
  rownames_to_column("micro")%>%
  mutate(micro =  factor(micro, levels = ordw12_micro)) %>%
  arrange(micro)%>%
  column_to_rownames("micro")
p<-pheatmap(corr_n, annotation_row = m16s_key_ext,annotation_col = mtb_key_ext, cluster_rows=FALSE, cluster_cols=FALSE,
            annotation_colors = ann_colors,color= inferno(10), main="Wk8")
ggsave(paste0(jrpca_res_wk8,"correlation_table/SFR24_0214_heatmap_wbdmannot_noclust.pdf"), p,width = 23, height = 25.5)
