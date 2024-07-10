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
library(EnhancedVolcano)
library(clusterProfiler)

#Transcriptomics Reanalysis of STAM-HCC data (ileum)

##################################################################

#clean metadata to just have NASH (wk12) samples

nash_score<-fread("transcriptomics/md_with_nash_score.tsv")%>%
  dplyr::select(sample_id,NASH_category)
md_il<-fread("transcriptomics/metadata.stam.noblanks.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  filter(sample_type=="ileum" & disease_state=="nash")%>%
  left_join(.,nash_score,by="sample_id")%>%
  dplyr::rename(sampleid=sample_id)%>%
  mutate(NASH_category=ifelse(condition=="NA","not_applicable",NASH_category))

write.table(md_il,"transcriptomics/ileum_md_with_nash_score.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

grade <- fread("files_from_KF/STAM_TRF/Histology/nas_scoring_analysis_nash_jingjing.csv", header=TRUE)%>%
  dplyr::select(mouse_id,steatosis_grade,fibrosis_stage)%>%
  dplyr::rename(host_subject_id=mouse_id)%>%
  mutate(steatosis_grade_new=ifelse(steatosis_grade>0, "1+","0"),
         fibrosis_stage_new=ifelse(fibrosis_stage>0, "1+","0"))

md_il<-md_il%>%left_join(.,grade, by="host_subject_id")%>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  mutate(steatosis_grade_new=ifelse(condition=="NA","none",steatosis_grade_new),
         fibrosis_stage_new=ifelse(condition=="NA","none",fibrosis_stage_new))

write.table(md_il,"transcriptomics/ileum_md_with_nash_score_fibste.tsv",sep = "\t",row.names = FALSE, quote=FALSE)


md_il_ZT1<-md_il%>%filter(timepoint=="ZT1")
md_il_ZT13<-md_il%>%filter(timepoint=="ZT13")
md_il_NA<-md_il%>%filter(condition=="NA")
md_il_FA<-md_il%>%filter(condition=="FA")
md_il_FT<-md_il%>%filter(condition=="FT")

#convert data file

dat<-fread("transcriptomics/data_files/ileum.nash.stringtie.gene_count_matrix.csv")
write.table(dat,"transcriptomics/data_files/ileum.nash.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#filter for just protein coding sequence
annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")%>%
  filter(molecule_type=="protein_coding")
dat_sub<-dat%>%filter(gene_id %in% annot$stringtie_id)
write.table(dat_sub,"transcriptomics/data_files/ileum.nash.stringtie.gene_count_prot_matrix.csv",sep = ",",row.names = FALSE, quote=FALSE)
write.table(dat_sub,"transcriptomics/data_files/ileum.nash.stringtie.gene_count_prot_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#filter for just FAFT protein coding sequence
md_il_FAFT<-md_il%>%filter(condition!="NA")
dat_FAFT<-dat_sub%>%dplyr::select(gene_id,all_of(md_il_FAFT$sampleid))
write.table(dat_FAFT,"transcriptomics/data_files/ileum.nash.stringtie.gene_count_FAFTprot_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat_FAFT,"transcriptomics/data_files/ileum.nash.stringtie.gene_count_FAFTprot_matrix.csv",sep = ",",row.names = FALSE, quote=FALSE)

#filter for just FAFT protein coding sequence & remove outlier 18L_I_S111
md_il_FAFT<-md_il%>%filter(condition!="NA" & sampleid!="18L_I_S111")
dat_FAFT<-dat_sub%>%dplyr::select(gene_id,all_of(md_il_FAFT$sampleid))
write.table(dat_FAFT,"transcriptomics/data_files/ileum.nash.stringtie.gene_count_FAFTprot_rmoutl_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat_FAFT,"transcriptomics/data_files/ileum.nash.stringtie.gene_count_FAFTprot_rmoutl_matrix.csv",sep = ",",row.names = FALSE, quote=FALSE)


#filter for just FA & remove outlier
md_il_FA<-md_il_FAFT%>%filter(condition=="FA")
dat_FA<-dat_FAFT%>%dplyr::select(gene_id,all_of(md_il_FA$sampleid))
write.table(dat_FA,"transcriptomics/data_files/ileum.nash.stringtie.gene_count_FAprot_rmoutl_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat_FA,"transcriptomics/data_files/ileum.nash.stringtie.gene_count_FAprot_rmoutl_matrix.csv",sep = ",",row.names = FALSE, quote=FALSE)

#filter for just FT
md_il_FT<-md_il_FAFT%>%filter(condition=="FT")
dat_FT<-dat_FAFT%>%dplyr::select(gene_id,all_of(md_il_FT$sampleid))
write.table(dat_FT,"transcriptomics/data_files/ileum.nash.stringtie.gene_count_FTprot_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat_FT,"transcriptomics/data_files/ileum.nash.stringtie.gene_count_FTprot_matrix.csv",sep = ",",row.names = FALSE, quote=FALSE)

dat_ZT1<-dat%>%dplyr::select(gene_id,md_il_ZT1$sampleid)
write.table(dat_ZT1,"transcriptomics/data_files/ileum.nash.ZT1.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_ZT13<-dat%>%dplyr::select(gene_id,md_il_ZT13$sampleid)
write.table(dat_ZT13,"transcriptomics/data_files/ileum.nash.ZT13.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_NA<-dat%>%dplyr::select(gene_id,md_il_NA$sampleid)
write.table(dat_NA,"transcriptomics/data_files/ileum.nash.NA.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_FA<-dat%>%dplyr::select(gene_id,md_il_FA$sampleid)
write.table(dat_FA,"transcriptomics/data_files/ileum.nash.FA.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_FT<-dat%>%dplyr::select(gene_id,md_il_FT$sampleid)
write.table(dat_FT,"transcriptomics/data_files/ileum.nash.FT.stringtie.gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

##################################################################

#NASH all
ord <- read_qza("transcriptomics/rpca_results_ileum/rpca_results_NASH_all/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"transcriptomics/rpca_results_ileum/rpca_results_NASH_all/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"transcriptomics/rpca_results_ileum/rpca_results_NASH_all/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("transcriptomics/ileum_md_with_nash_score_fibste.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(sample_name=sampleid)

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(timepoint=="ZT1","light","dark"))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","not_applicable")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=NASH_category)) +
  geom_point(alpha=1.0, size=2.5) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = condition_ZT))+
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_shape_manual(values=c(16,17,3)) +
  #scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH ileum RNA (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("transcriptomics/rpca_results_ileum/rpca_results_NASH_all/SFR23_1024_NASH_ileumRNA_RPCA.pdf", plot=p,height=4, width=4)
ggsave("transcriptomics/rpca_results_ileum/rpca_results_NASH_all/SFR23_1024_NASH_ileumRNA_NASHcat_RPCA.pdf", plot=p,height=4, width=4)

#NASH all (prot)
ord <- read_qza("transcriptomics/rpca_results_ileum/rpca_results_NASH_prot/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"transcriptomics/rpca_results_ileum/rpca_results_NASH_prot/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"transcriptomics/rpca_results_ileum/rpca_results_NASH_prot/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("transcriptomics/ileum_md_with_nash_score_fibste.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(sample_name=sampleid)

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(timepoint=="ZT1","light","dark"))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","not_applicable")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=NASH_category)) +
  geom_point(alpha=1.0, size=2.5) + #alpha controls transparency and helps when points are overlapping
  theme_classic() +#stat_ellipse(type = "t", linetype = 2,aes(group = condition_ZT))+
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_shape_manual(values=c(16,17,3)) +
  #scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH ileum RNA (Wk 12)")+ 
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")
ggsave("transcriptomics/rpca_results_ileum/rpca_results_NASH_prot/SFR23_1212_NASH_ileumRNA_NASHcat_RPCA.pdf", plot=p,height=4, width=4)

#NASH FAFT (prot)
ord <- read_qza("transcriptomics/rpca_results_ileum/rpca_results_NASH_FAFTprot/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"transcriptomics/rpca_results_ileum/rpca_results_NASH_FAFTprot/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"transcriptomics/rpca_results_ileum/rpca_results_NASH_FAFTprot/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("transcriptomics/ileum_md_with_nash_score_fibste.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(sample_name=sampleid)

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(timepoint=="ZT1","light","dark"))%>%
  mutate(condition=factor(condition,levels=c("FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=NASH_category)) +
  geom_point(alpha=1.0, size=2.5) + #alpha controls transparency and helps when points are overlapping
  theme_classic() +#stat_ellipse(type = "t", linetype = 2,aes(group = condition_ZT))+
  scale_color_manual(values=c("#D55E00","#009E73"))+
  scale_shape_manual(values=c(16,17)) +
  #scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH ileum RNA (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("transcriptomics/rpca_results_ileum/rpca_results_NASH_FAFTprot/SFR23_1212_NASH_ileumRNA_NASHcat_RPCA.pdf", plot=p,height=4, width=4)

##################################################################

get_genelist<-function(df){
  annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")%>%
    dplyr::select(stringtie_id,ncbi_id)
  df<-df%>%dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")
  
  original_gene_list <- df$log2FoldChange
  names(original_gene_list) <- df$ncbi_id
  gene_list<-na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  return(gene_list)
}

gse_nash_plot<-function(df){
  set.seed(1234)
  gene_list<-get_genelist(df)
  # gse <- gseGO(geneList=gene_list,
  #              ont ="ALL",
  #              keyType = "ENTREZID",
  #              minGSSize = 3,
  #              maxGSSize = 800,
  #              pvalueCutoff = 0.05,
  #              verbose = TRUE,
  #              OrgDb = organism,
  #              pAdjustMethod = "none")
  
  kk2 <- gseKEGG(geneList     = gene_list,
                 organism     = "mmu",
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 keyType       = "ncbi-geneid")

  #p<-dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  p<-dotplot(kk2, showCategory=10, split=".sign") + facet_grid(.~.sign)
  return(p)
}

gse_compareclust<-function(df1,df2,cond1,cond2){
  
  gene_list1<-get_genelist(df1)
  gene_list2<-get_genelist(df2)
  
  mydf1 <- data.frame(Entrez=names(gene_list1), FC=gene_list1)%>%
    mutate(direction=ifelse(FC<0,"a-Non-NASH","b-NASH"),
           #mutate(direction=ifelse(FC<0,"FA","FT"),
           #mutate(direction=ifelse(FC<0,"ZT1","ZT13"),
           condition=cond1)
  
  mydf2 <- data.frame(Entrez=names(gene_list2), FC=gene_list2)%>%
    mutate(direction=ifelse(FC<0,"a-Non-NASH","b-NASH"),
           #mutate(direction=ifelse(FC<0,"FA","FT"),
           #mutate(direction=ifelse(FC<0,"ZT1","ZT13"),
           condition=cond2)
  
  mydf<-rbind(mydf1,mydf2)
  # formula_res <- compareCluster(Entrez~direction+condition, OrgDb = "org.Mm.eg.db", data=mydf,
  #                               fun="enrichGO",pvalueCutoff=0.05,pAdjustMethod="fdr")
  formula_res <- compareCluster(Entrez~direction+condition, organism = "mmu", data=mydf,
                                fun="enrichKEGG",pvalueCutoff=0.05,pAdjustMethod="fdr")
  return(formula_res)
}

gse_compareclust_1<-function(df,cond){
  
  gene_list<-get_genelist(df)
  
  mydf <- data.frame(Entrez=names(gene_list), FC=gene_list)%>%
    mutate(direction=ifelse(FC<0,"a-Non-NASH","b-NASH"),
           condition=cond)
  
  # formula_res <- compareCluster(Entrez~direction+condition, OrgDb = "org.Mm.eg.db", data=mydf,
  #                               fun="enrichGO",pvalueCutoff=0.05,pAdjustMethod="fdr")
  formula_res <- compareCluster(Entrez~direction+condition, organism = "mmu", data=mydf,
                                fun="enrichKEGG",pvalueCutoff=0.05,pAdjustMethod="fdr")
  return(formula_res)
}

#nash_category

organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
require(DOSE)

#log(NASH/Non-NASH)
deq_nash<-fread("transcriptomics/NASH_ileum_deseq2/FAFTprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-NASH-Non_NASH-pval_05.csv")
gse_nash_plot(deq_nash)+ggtitle("Non-NASH vs NASH")+scale_y_discrete(position = "right")
#ggsave("transcriptomics/NASH_ileum_deseq2/FAFTprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-NASH-Non_NASH.gse.pdf",height=5, width=8)
ggsave("transcriptomics/NASH_ileum_deseq2/FAFTprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-NASH-Non_NASH.kegg.pdf",height=4, width=9)

il_nash<-gse_compareclust_1(deq_nash,"both")
il_nash_df<-il_nash@compareClusterResult%>%
  mutate(direction=factor(direction,levels=c("a-Non-NASH","b-NASH")))%>%
  filter(p.adjust<0.05)
il_nash@compareClusterResult = il_nash@compareClusterResult[il_nash@compareClusterResult$ID %in% il_nash_df$ID, ]

dotplot(il_nash,x="direction", showCategory=15) + 
  facet_grid(category~., scales="free_y",space='free_y') + scale_y_discrete(position = "right") + 
  theme(axis.text.y = element_text(size=8),
        legend.position = "left") +ggtitle("Ileum")+
  scale_x_discrete(labels= c("a-Non-NASH","b-NASH"))
ggsave("transcriptomics/NASH_ileum_deseq2/FAFTprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-NASH-Non_NASH.kegg.comparecluster.pdf",height=6, width=7)

# deq_nash_time<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.timepoint.NASH_category-ZT1.NASH-ZT13.NASH-pval_05.csv")
# gse_nash_plot(deq_nash_time)+ggtitle("ZT1 NASH vs ZT13 NASH")
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.NASH_category-ZT1.NASH-ZT13.NASH.gse.pdf",height=8, width=7)
# 
# deq_nash_time_ZT13<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.timepoint.NASH_category-ZT13.NASH-ZT13.Non_NASH-pval_05.csv")
# gse_nash_plot(deq_nash_time_ZT13)+ggtitle("Non-NASH vs NASH ZT13")
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.NASH_category-ZT13.NASH-ZT13.Non_NASH.gse.pdf",height=8, width=7)


deq_nash_FA<-fread('transcriptomics/NASH_ileum_deseq2/FAprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-FA.NASH-FA.Non_NASH-pval_05.csv')
gse_nash_plot(deq_nash_FA)+ggtitle("FA non-NASH vs NASH")+scale_y_discrete(position = "right")
#ggsave("transcriptomics/NASH_ileum_deseq2/FAprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-FA.NASH-FA.Non_NASH.gse.pdf",height=5, width=10)
ggsave("transcriptomics/NASH_ileum_deseq2/FAprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-FA.NASH-FA.Non_NASH.kegg.pdf",height=5, width=10)

deq_nash_FT<-fread('transcriptomics/NASH_ileum_deseq2/FTprot/NASH_ileum_deseq2_.lfchange.NASH_category-FT.NASH-FT.Non_NASH-pval_05.csv')
gse_nash_plot(deq_nash_FT)+ggtitle("FT non-NASH vs NASH")+scale_y_discrete(position = "right")
#ggsave("transcriptomics/NASH_ileum_deseq2/FTprot/NASH_ileum_deseq2_.lfchange.NASH_category-FT.NASH-FT.Non_NASH.gse.pdf",height=5, width=8)
ggsave("transcriptomics/NASH_ileum_deseq2/FTprot/NASH_ileum_deseq2_.lfchange.NASH_category-FT.NASH-FT.Non_NASH.kegg.pdf",height=2, width=8)

#plot together using comparecluster
deq_nash_FA<-fread('transcriptomics/NASH_ileum_deseq2/FAprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-FA.NASH-FA.Non_NASH-pval_05.csv')
deq_nash_FT<-fread('transcriptomics/NASH_ileum_deseq2/FTprot/NASH_ileum_deseq2_.lfchange.NASH_category-FT.NASH-FT.Non_NASH-pval_05.csv')

FAFT_nash<-gse_compareclust(deq_nash_FA,deq_nash_FT,"FA","FT")

dotplot(FAFT_nash,x="direction", showCategory=15) + 
  facet_grid(~condition) + scale_y_discrete(position = "right") + 
  theme(axis.text.y = element_text(size=8),
        legend.position = "left") +ggtitle("Ileum")
ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.NASH_category-nashnonnash-FAFT.compclust.pdf",height=6.5, width=8)

dotplot(FAFT_nash,x="direction", showCategory=50) + 
  facet_grid(category~condition, scales="free_y",space='free_y') + scale_y_discrete(position = "right") + 
  theme(axis.text.y = element_text(size=8),
        legend.position = "left") +ggtitle("Ileum")
ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.NASH_category-nashnonnash-FAFT.compclust.kegg.pdf",height=4, width=8)

FAFT_nash_df<-FAFT_nash%>%as.data.frame()

#get all the genes with transmembrane activity under non-nash
annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")

# NASH_FA_chn<-FAFT_nash_df%>%filter(Cluster=="a-Non-NASH.FA" & grepl("channel",Description))%>%
#   separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
#   dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
#   left_join(.,annot,by="ncbi_id")%>%
#   dplyr::rename(gene_id=stringtie_id)

nonNASH_ppar<-il_nash_df%>%filter(direction=="a-Non-NASH" & grepl("PPAR signaling pathway",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="PPAR")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)


NASH_FA_trans<-FAFT_nash_df%>%filter(Cluster=="a-Non-NASH.FA" & grepl("transmembrane transporter",Description))%>% #(carboxylic acid|amino acid|organic acid) 
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  dplyr::rename(gene_id=stringtie_id)

md<-fread("transcriptomics/ileum_md_with_nash_score_fibste.tsv") %>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))


#plot heatmap of these specific hits

#genecounts<-fread("transcriptomics/data_files/stringtie.gene_tpm_matrix.csv")%>%
genecounts<-fread("transcriptomics/NASH_ileum_deseq2/FAFTprot_rmoutlr/NASH_ileum_deseq2_.vsd.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  dplyr::select(gene_id,any_of(md$sampleid))%>%
  #filter(gene_id %in% NASH_FA_chn$gene_id)%>%
  #filter(gene_id %in% NASH_FA_trans$gene_id)%>%
  filter(gene_id %in% nonNASH_ppar$gene_id)%>%
  gather(sampleid,TPM_counts,-gene_id) %>%
  left_join(.,md,by="sampleid") %>%
  filter(condition!="NA")%>%
  #left_join(.,NASH_FA_chn,by="gene_id")%>%
  #left_join(.,NASH_FA_trans,by="gene_id")%>%
  left_join(.,nonNASH_ppar,by="gene_id")%>%
  group_by(gene_id,condition,NASH_category)%>%dplyr::summarise(mn_TPM=mean(TPM_counts))%>%
  group_by(gene_id)%>%mutate(Zscore=(mn_TPM - mean(mn_TPM))/sd(mn_TPM))%>%
  mutate(condition=factor(condition,levels=c("FA","FT")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","not_applicable")))%>%
  filter(!is.na(NASH_category))%>%
  separate(gene_id, c("gene_number","gene_symbol"), remove=FALSE)


# Run clustering
genecounts_df<-fread("transcriptomics/NASH_ileum_deseq2/FAFTprot_rmoutlr/NASH_ileum_deseq2_.vsd.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  dplyr::select(gene_id,any_of(md$sampleid))%>%
  #filter(gene_id %in% NASH_FA_trans$gene_id)
  filter(gene_id %in% nonNASH_ppar$gene_id)
genecounts_matrix<-genecounts_df%>%column_to_rownames("gene_id") %>%as.matrix()
genecounts_dendro <- as.dendrogram(hclust(d = dist(x = genecounts_matrix)))
dendro_plot <- ggdendrogram(data = genecounts_dendro, rotate = TRUE)

genecounts_order <- order.dendrogram(genecounts_dendro)
genecounts$gene_id <- factor(x = genecounts$gene_id,
                             levels = genecounts_df$gene_id[genecounts_order], 
                             ordered = TRUE)


plt<-ggplot(genecounts,aes(x=NASH_category, y=gene_symbol)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(~condition,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),axis.text.y = element_text(size = 7),
        panel.spacing.y=unit(0.01, "lines"),strip.text.y = element_text(angle = 0))+
  scale_fill_viridis_b(option = "plasma")

#ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_channel_hits_GOenrich_v2.pdf", plot=plt,height=6, width=4)
#ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_trans_hits_GOenrich_v2.pdf", plot=plt,height=4, width=4)
ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_ppar_hits_KEGGenrich.pdf", plot=plt,height=2, width=4)


# cnetplot(FAFT_nash)+scale_fill_manual(values=c("#E69F00","limegreen","#D55E00","darkgreen"))
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.NASH_category-nashnonnash-FAFT-network.compclust.pdf",height=15, width=15)
# 
# deq_nash_nash<-fread('transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.NASH_category-FA.NASH-FT.NASH-pval_05.csv')
# gse_nash_plot(deq_nash_nash)+ggtitle("FT vs. FA NASH")
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.NASH_category-FA.NASH-FT.NASH.gse.pdf",height=8.5, width=8)

# deq_nash_nonnash<-fread('transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.NASH_category-FA.Non_NASH-FT.Non_NASH-pval_05.csv')
# gse_nash_plot(deq_nash_nonnash)+ggtitle("FA vs. FT non-NASH")
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.NASH_category-FA.Non_NASH-FT.Non_NASH.gse.pdf",height=8, width=7)
# 
# deq_nash_nash<-fread('transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.NASH_category-FA.NASH-FT.NASH-pval_05.csv')
# deq_nash_nonnash<-fread('transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.NASH_category-FA.Non_NASH-FT.Non_NASH-pval_05.csv')
# nonnash_nash<-gse_compareclust(deq_nash_nonnash,deq_nash_nash,"a-Non-NASH","b-NASH")

# dotplot(nonnash_nash,x="direction", showCategory=15) + facet_grid(~condition)
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.NASH_category-FAFT-nashnonnash.compclust.pdf",height=10, width=10)
# 
# cnetplot(nonnash_nash)+scale_fill_manual(values=c("#E69F00","#D55E00","limegreen","darkgreen"))
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.NASH_category-FAFT-nashnonnash-network.compclust.pdf",height=15, width=15)

# deq_timept<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-NA.ZT1-NA.ZT13-pval_05.csv")
# gse_nash_plot(deq_timept)+ggtitle("ZT1 vs ZT13 NA")
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.timepoint-NA.ZT1-NA.ZT13.gse.pdf",height=8, width=7)
# 
# deq_timept_FA<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-FA.ZT1-FA.ZT13-pval_05.csv")
# gse_nash_plot(deq_timept_FA)+ggtitle("ZT1 vs ZT13 FA")
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.timepoint-FA.ZT1-FA.ZT13.gse.pdf",height=8, width=7)
# 
# deq_timept_FT<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-FT.ZT1-FT.ZT13-pval_05.csv")
# gse_nash_plot(deq_timept_FT)+ggtitle("ZT1 vs ZT13 FT")
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.timepoint-FT.ZT1-FT.ZT13.gse.pdf",height=8, width=7)
# 
# deq_timept_FA<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-FA.ZT1-FA.ZT13-pval_05.csv")
# deq_timept_FT<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-FT.ZT1-FT.ZT13-pval_05.csv")
# ZT_FAFT<-gse_compareclust(deq_timept_FA,deq_timept_FA,"FA","FT")
# 
# dotplot(ZT_FAFT,x="direction", showCategory=25) + facet_grid(~condition)
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.NASH_category-FT.ZT1-FT.ZT13.compclust.pdf",height=10, width=10)
# 
# cnetplot(ZT_FAFT)+scale_fill_manual(values=c("#E69F00","limegreen","#D55E00","darkgreen"))
# ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2.NASH_category-FT.ZT1-FT.ZT13-network.compclust.pdf",height=15, width=15)


deq_nash_FA<-fread('transcriptomics/NASH_ileum_deseq2/FAprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-FA.NASH-FA.Non_NASH-pval_05.csv')
deq_nash_FT<-fread('transcriptomics/NASH_ileum_deseq2/FTprot/NASH_ileum_deseq2_.lfchange.NASH_category-FT.NASH-FT.Non_NASH-pval_05.csv')
list_venn <- list(FA = deq_nash_FA$V1,
                  FT = deq_nash_FT$V1)

p<-venn.diagram(list_venn, fill = c("#D55E00","#009E73"),height = 10,lty = 0, 
                main="non-NASH vs NASH",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_venn_nonnash_nash.pdf")
#pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_venn_non-nash.pdf")
#pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_venn_nash.pdf")
grid.draw(p)
dev.off()

# deq_nash_nash<-fread('transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.NASH_category-FA.NASH-FT.NASH-pval_05.csv')
# deq_nash_nonnash<-fread('transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.NASH_category-FA.Non_NASH-FT.Non_NASH-pval_05.csv')
# 
# list_venn <- list(NASH = deq_nash_nash$V1,
#                   NonNASH = deq_nash_nonnash$V1)
# 
# p<-venn.diagram(list_venn, fill = c("#D55E00","#009E73"),height = 10,lty = 0, 
#                 main="FA vs. FT",
#                 width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)
# 
# grid.draw(p)
# pdf(file="transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_venn_FAFT_nonnash_nash.pdf")
# grid.draw(p)
# dev.off()

# deq_timept_FA<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-FA.ZT1-FA.ZT13-pval_05.csv")
# deq_timept_FT<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-FT.ZT1-FT.ZT13-pval_05.csv")
# 
# list_venn <- list(FA = deq_timept_FA$V1,
#                   FT = deq_timept_FT$V1)
# 
# p<-venn.diagram(list_venn, fill = c("#D55E00","#009E73"),height = 10,lty = 0, 
#                 main="ZT1 vs ZT13",
#                 width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)
# 
# grid.draw(p)
# pdf(file="transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_venn_ZT1ZT13.pdf")
# grid.draw(p)
# dev.off()

##################################################################

#draw histogram to see distribution
# histo_dat<-fread("transcriptomics/data_files/ileum.nash.stringtie.gene_count_prot_matrix.csv")%>%
#   gather(sample_name,counts,-gene_id)
# 
# ggplot(histo_dat, aes(x=counts)) + geom_density() + coord_cartesian(xlim=c(-5,2500))

annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")%>%
  filter(molecule_type=="protein_coding")

deq_nash<-fread("transcriptomics/NASH_ileum_deseq2/FAFTprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-NASH-Non_NASH.csv")%>%
  dplyr::rename(stringtie_id=V1)%>%
  left_join(.,annot,by="stringtie_id")%>%
  filter(!is.na(padj))

p<-EnhancedVolcano(
  deq_nash,
  lab = deq_nash$gene_symbol,
  title = 'non-NASH vs. NASH',
  x = "log2FoldChange",
  y = "padj",
  labSize = 3.0,
  subtitle=NA,
  caption=NA,
  xlim=c(-4,4),
  ylim=c(0,10),
  colAlpha=1,
  FCcutoff = 1,
  pointSize = 1,
  legendPosition = "none",
  pCutoff = 0.001,
  gridlines.major=FALSE,
  gridlines.minor=FALSE,
  border="full",
) +theme(text=element_text(size=7))

ggsave("transcriptomics/NASH_ileum_deseq2/FAFTprot_rmoutlr/NASH_ileum_deseq2.NASH_category-NASH-Non_NASH.volplot.pdf",height=5, width=5.5)

annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")%>%
  filter(molecule_type=="protein_coding")

deq_nash_FA<-fread('transcriptomics/NASH_ileum_deseq2/FAprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-FA.NASH-FA.Non_NASH.csv')%>%
  mutate(log2FoldChange=log2FoldChange)%>%
  dplyr::rename(stringtie_id=V1)%>%
  left_join(.,annot,by="stringtie_id")%>%
  filter(!is.na(padj))

p<-EnhancedVolcano(
  deq_nash_FA,
  lab = deq_nash_FA$gene_symbol,
  title = 'non-NASH vs. NASH (FA)',
  x = "log2FoldChange",
  y = "padj",
  labSize = 3.0,
  subtitle=NA,
  caption=NA,
  colAlpha=1,
  FCcutoff = 1,
  pointSize = 1,
  col=c('gray37', 'gray37','darkgoldenrod2','darkorange3'),
  legendPosition = "none",
  pCutoff = 0.001,
  gridlines.major=FALSE,
  gridlines.minor=FALSE,
  border="full",
) +theme(text=element_text(size=7))

ggsave("transcriptomics/NASH_ileum_deseq2/FAprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-FA.NASH-FA.Non_NASH.volplot.pdf",height=5, width=5)


deq_nash_FT<-fread('transcriptomics/NASH_ileum_deseq2/FTprot/NASH_ileum_deseq2_.lfchange.NASH_category-FT.NASH-FT.Non_NASH.csv')%>%
  mutate(log2FoldChange=-1*log2FoldChange)%>%
  dplyr::rename(stringtie_id=V1)%>%
  left_join(.,annot,by="stringtie_id")%>%
  filter(!is.na(padj))

p<-EnhancedVolcano(
  deq_nash_FT,
  lab = deq_nash_FT$gene_symbol,
  title = 'non-NASH vs. NASH (FT)',
  x = "log2FoldChange",
  y = "padj",
  labSize = 3.0,
  subtitle=NA,
  caption=NA,
  colAlpha=1,
  FCcutoff = 1,
  pointSize = 1,
  col=c('gray37', 'gray37', "darkgreen","chartreuse3"),
  legendPosition = "none",
  pCutoff = 0.001,
  gridlines.major=FALSE,
  gridlines.minor=FALSE,
  border="full",
) +theme(text=element_text(size=7))

ggsave("transcriptomics/NASH_ileum_deseq2/FTprot/NASH_ileum_deseq2_.lfchange.NASH_category-FT.NASH-FT.Non_NASH.volplot.pdf",height=5, width=5)

##################################################################

#make heatmap of functions of interest
clock_genes<-c("Arntl","Clock","Cry2","Dbp","Nr1d1","Nr1d2","Per1","Per2","Rora")
BA_genes<-c("Fabp6","Slc10a2","Fgf15","Nr0b2","Nr1h4","Slc51b")
glp_genes<-c("Dpp4","Ffar2","Ffar4","Gcg","Slc2a2","Slc2a5")
annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")

deq_nash_FA<-fread('transcriptomics/NASH_ileum_deseq2/FAprot_rmoutlr/NASH_ileum_deseq2_.lfchange.NASH_category-FA.NASH-FA.Non_NASH.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")
deq_nash_FT<-fread('transcriptomics/NASH_ileum_deseq2/FTprot/NASH_ileum_deseq2_.lfchange.NASH_category-FT.NASH-FT.Non_NASH.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")

clk_FA<-deq_nash_FA%>%filter(gene_symbol %in% clock_genes)%>%mutate(condition="FA")
clk_FT<-deq_nash_FT%>%filter(gene_symbol %in% clock_genes)%>%mutate(condition="FT")

clk<-rbind(clk_FA,clk_FT)

plt<-ggplot(clk,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(color="black",aes(fill=log2FoldChange))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000", #colors in the scale
                       midpoint=0,    #same midpoint for plots (mean of the range)
                       breaks=seq(-100,100,2), #breaks in the scale bar
                       limits=c(-5,5))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  geom_text(aes(label = signif(padj,digits=2)), color = "black") +
  guides(fill = guide_colourbar(barheight = 0.5))+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0), legend.positio="top") + 
  labs(x="",y="")+ggtitle("Circadian Clock Genes\nNon-NASH v. NASH") 

ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.NASH_category-clockgenes.pdf", plot=plt,height=4, width=3)


BA_FA<-deq_nash_FA%>%filter(gene_symbol %in% BA_genes)%>%mutate(condition="FA")
BA_FT<-deq_nash_FT%>%filter(gene_symbol %in% BA_genes)%>%mutate(condition="FT")

BA<-rbind(BA_FA,BA_FT)

plt<-ggplot(BA,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(color="black", aes(fill=log2FoldChange))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000", #colors in the scale
                       midpoint=0,    #same midpoint for plots (mean of the range)
                       breaks=seq(-100,100,1), #breaks in the scale bar
                       limits=c(-3,3))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  geom_text(aes(label = signif(padj,digits=2)), color = "black") +
  guides(fill = guide_colourbar(barheight = 0.5))+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0),legend.position = "top") + 
  labs(x="",y="")+ggtitle("Bile Acid Genes\nNon-NASH v. NASH")

ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.NASH_category-BAgenes.pdf", plot=plt,height=4, width=3)


glp_FA<-deq_nash_FA%>%filter(gene_symbol %in% glp_genes)%>%mutate(condition="FA")
glp_FT<-deq_nash_FT%>%filter(gene_symbol %in% glp_genes)%>%mutate(condition="FT")

glp<-rbind(glp_FA,glp_FT)

plt<-ggplot(BA,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(color="black", aes(fill=log2FoldChange))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000", #colors in the scale
                       midpoint=0,    #same midpoint for plots (mean of the range)
                       breaks=seq(-100,100,1), #breaks in the scale bar
                       limits=c(-3,3))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  geom_text(aes(label = signif(padj,digits=2)), color = "black") +
  guides(fill = guide_colourbar(barheight = 0.5))+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0),legend.position = "top") + 
  labs(x="",y="")+ggtitle("GLP-1 Signaling\nNon-NASH v. NASH")

ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.NASH_category-glpgenes.pdf", plot=plt,height=4, width=3)


deq_timept_FA<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-FA.ZT1-FA.ZT13.csv")%>%
  dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")
deq_timept_FT<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-FT.ZT1-FT.ZT13.csv")%>%
  dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")
deq_timept_NA<-fread("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-NA.ZT1-NA.ZT13.csv")%>%
  dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")

clk_FA<-deq_timept_FA%>%filter(gene_symbol %in% clock_genes)%>%mutate(condition="FA")
clk_FT<-deq_timept_FT%>%filter(gene_symbol %in% clock_genes)%>%mutate(condition="FT")
clk_NA<-deq_timept_NA%>%filter(gene_symbol %in% clock_genes)%>%mutate(condition="NA")

clk<-rbind(clk_FA,clk_FT,clk_NA)%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")))

plt<-ggplot(clk,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(aes(fill=log2FoldChange))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0),legend.position = "top") + 
  scale_fill_viridis_c() + labs(x="",y="")+ggtitle("Circadian Clock Genes\nZT1 vs. ZT13")

ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-clockgenes.pdf", plot=plt,height=3.5, width=2.5)

BA_FA<-deq_timept_FA%>%filter(gene_symbol %in% BA_genes)%>%mutate(condition="FA")
BA_FT<-deq_timept_FT%>%filter(gene_symbol %in% BA_genes)%>%mutate(condition="FT")
BA_NA<-deq_timept_NA%>%filter(gene_symbol %in% BA_genes)%>%mutate(condition="NA")

BA<-rbind(BA_FA,BA_FT,BA_NA)%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")))

plt<-ggplot(BA,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(aes(fill=log2FoldChange))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0),legend.position = "top") + 
  scale_fill_viridis_c() + labs(x="",y="")+ggtitle("Bile Acid Genes\nZT1 vs. ZT13")

ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-BAgenes.pdf", plot=plt,height=3.5, width=2.5)

glp_FA<-deq_timept_FA%>%filter(gene_symbol %in% glp_genes)%>%mutate(condition="FA")
glp_FT<-deq_timept_FT%>%filter(gene_symbol %in% glp_genes)%>%mutate(condition="FT")
glp_NA<-deq_timept_NA%>%filter(gene_symbol %in% glp_genes)%>%mutate(condition="NA")

glp<-rbind(glp_FA,glp_FT,glp_NA)%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")))

plt<-ggplot(glp,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(aes(fill=log2FoldChange))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0),legend.position = "top") + 
  scale_fill_viridis_c() + labs(x="",y="")+ggtitle("GLP-1 signaling\nZT1 vs. ZT13")

ggsave("transcriptomics/NASH_ileum_deseq2/NASH_ileum_deseq2_.lfchange.condition.timepoint-glpgenes.pdf", plot=plt,height=3.5, width=2.5)

