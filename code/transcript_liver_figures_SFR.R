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
library(clusterProfiler)
library(enrichplot)
library(EnhancedVolcano)
library(viridis)
library("ggdendro")
library(cowplot)

#Transcriptomics Reanalysis of STAM-HCC data (liver)
##################################################################

md_liv<-fread("transcriptomics/liver_metadata.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  dplyr::rename(sampleid=sample_id,host_subject_id=host_id)

grade <- fread("files_from_KF/STAM_TRF/Histology/nas_scoring_analysis_nash_jingjing.csv", header=TRUE)%>%
  dplyr::select(mouse_id,steatosis_grade,fibrosis_stage,NASH_category)%>%
  dplyr::rename(host_subject_id=mouse_id)%>%
  mutate(steatosis_grade_new=ifelse(steatosis_grade>0, "1+","0"),
         fibrosis_stage_new=ifelse(fibrosis_stage>0, "1+","0"))

md_liv<-md_liv%>%left_join(.,grade, by="host_subject_id")%>%
  mutate(steatosis_grade_new=ifelse(condition=="NA","none",steatosis_grade_new),
         fibrosis_stage_new=ifelse(condition=="NA","none",fibrosis_stage_new),
         NASH_category=ifelse(condition=="NA","not_applicable",NASH_category))

write.table(md_liv,"transcriptomics/liver_metadata_wnashscore_fibste.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

md_liv_ZT1<-md_liv%>%filter(timepoint=="ZT1")
md_liv_ZT13<-md_liv%>%filter(timepoint=="ZT13")
md_liv_NA<-md_liv%>%filter(condition=="NA")
md_liv_FA<-md_liv%>%filter(condition=="FA")
md_liv_FT<-md_liv%>%filter(condition=="FT")

#convert data file
dat<-fread("transcriptomics/data_files/liver_gene_count_matrix.csv")
write.table(dat,"transcriptomics/data_files/liver_gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#filter for just protein coding sequence
annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")%>%
  filter(molecule_type=="protein_coding")
dat_sub<-dat%>%filter(gene_id %in% annot$stringtie_id)
write.table(dat_sub,"transcriptomics/data_files/liver_gene_count_prot_matrix.csv",sep = ",",row.names = FALSE, quote=FALSE)
write.table(dat_sub,"transcriptomics/data_files/liver_gene_count_prot_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#filter for just FA and FT (prot coding)
md_liv_FAFT<-md_liv%>%filter(condition!="NA")
dat_FAFT<-dat_sub%>%dplyr::select(gene_id,all_of(md_liv_FAFT$sampleid))
write.table(dat_FAFT,"transcriptomics/data_files/liver_gene_count_FAFTprot_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_ZT1<-dat%>%dplyr::select(gene_id,md_liv_ZT1$sampleid)
write.table(dat_ZT1,"transcriptomics/data_files/liver_ZT1_gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_ZT13<-dat%>%dplyr::select(gene_id,md_liv_ZT13$sampleid)
write.table(dat_ZT13,"transcriptomics/data_files/liver_ZT13_gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_NA<-dat%>%dplyr::select(gene_id,md_liv_NA$sampleid)
write.table(dat_NA,"transcriptomics/data_files/liver_NA_gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_FA<-dat%>%dplyr::select(gene_id,md_liv_FA$sampleid)
write.table(dat_FA,"transcriptomics/data_files/liver_FA_gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_FT<-dat%>%dplyr::select(gene_id,md_liv_FT$sampleid)
write.table(dat_FT,"transcriptomics/data_files/liver_FT_gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
##################################################################
#plot rpca

#NASH all
ord <- read_qza("transcriptomics/rpca_results_liver/rpca_results_NASH_all/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"transcriptomics/rpca_results_liver/rpca_results_NASH_all/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"transcriptomics/rpca_results_liver/rpca_results_NASH_all/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("transcriptomics/liver_metadata_wnashscore_fibste.tsv")%>%
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
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH liver RNA (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("transcriptomics/rpca_results_liver/rpca_results_NASH_all/SFR23_1024_NASH_liverRNA_RPCA.pdf", plot=p,height=4, width=4)
ggsave("transcriptomics/rpca_results_liver/rpca_results_NASH_all/SFR23_1024_NASH_liverRNA_NASHcat_RPCA.pdf", plot=p,height=4, width=4)

#NASH all (prot only)
ord <- read_qza("transcriptomics/rpca_results_liver/rpca_results_NASH_prot/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"transcriptomics/rpca_results_liver/rpca_results_NASH_prot/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"transcriptomics/rpca_results_liver/rpca_results_NASH_prot/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("transcriptomics/liver_metadata_wnashscore_fibste.tsv")%>%
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
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH liver RNA (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("transcriptomics/rpca_results_liver/rpca_results_NASH_prot/SFR23_1212_NASH_liverRNA_NASHcat_RPCA.pdf", plot=p,height=4, width=4)

#NASH FAFT (prot only)
ord <- read_qza("transcriptomics/rpca_results_liver/rpca_results_NASH_FAFTprot/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"transcriptomics/rpca_results_liver/rpca_results_NASH_FAFTprot/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"transcriptomics/rpca_results_liver/rpca_results_NASH_FAFTprot/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("transcriptomics/liver_metadata_wnashscore_fibste.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(sample_name=sampleid)%>%
  mutate(ZT_nash=paste(timepoint,NASH_category, sep="_"))%>%
  mutate(ZT_nash=factor(ZT_nash,levels=c("ZT1_NASH","ZT13_Non_NASH","ZT13_NASH")))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(timepoint=="ZT1","light","dark"))%>%
  mutate(condition=factor(condition,levels=c("FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, shape=NASH_category)) +
  geom_point(aes(fill=condition),alpha=1.0, size=2.5) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = condition_ZT))+
  scale_fill_manual(values=c("#D55E00","#009E73"))+
  scale_shape_manual(values=c(21,24)) +
  labs(x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH liver RNA (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("transcriptomics/rpca_results_liver/rpca_results_NASH_FAFTprot/SFR23_1212_NASH_liverRNA_NASHcat_RPCA_v2.pdf", plot=p,height=4, width=4)

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, shape=ZT_nash)) +
  geom_point(aes(fill=condition),alpha=1.0, size=2.5) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = condition_ZT))+
  scale_fill_manual(values=c("#D55E00","#009E73"))+
  scale_shape_manual(values=c(21,22,24)) +
  labs(x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH liver RNA (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("transcriptomics/rpca_results_liver/rpca_results_NASH_FAFTprot/SFR24_0417_NASH_liverRNA_ZTnash_RPCA_v2.pdf", plot=p,height=4, width=4)

##add box plot to axes

cond_rpca <- ggplot(rpca, aes(x =condition, y = PC2,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_point(colour="black", aes(fill=condition),
                                       position=position_jitter(width=0.1,height = 0.1),
                                       size=3,pch=21) + 
  labs(x = NA, y = "PC2") +
  scale_fill_manual(values=c("#D55E00","#009E73"))+
  scale_color_manual(values=c("#D55E00","#009E73"))+
  # theme_void()
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

pairwise.t.test(rpca$PC2, rpca$condition, p.adjust.method = "fdr")
# FA     
# FT 9.7e-07

nash_rpca <- ggplot(rpca, aes(x =NASH_category, y = PC1,shape=NASH_category)) +
  geom_boxplot(alpha=0.3) + geom_point(fill="black",colour="white",position=position_jitter(width=0.1,height = 0.1),size=3) +
  scale_shape_manual(values=c(21,24)) +
  scale_x_discrete(limits=rev)+
  labs(x = NA, y = "PC1") +
  # theme_void()
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+ coord_flip()

pairwise.t.test(rpca$PC1, rpca$NASH_category, p.adjust.method = "fdr")
# Non_NASH
# NASH 1.6e-05 

ZT_rpca <- ggplot(rpca, aes(x =ZT_nash, y = PC1,shape=ZT_nash)) +
  geom_boxplot(alpha=0.3) + geom_point(fill="black",colour="white",position=position_jitter(width=0.1,height = 0.1),size=3) +
  scale_shape_manual(values=c(21,22,24)) +
  scale_x_discrete(limits=rev)+
  labs(x = NA, y = "PC1") +
  # theme_void()
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+ coord_flip()

pairwise.t.test(rpca$PC1, rpca$timepoint, p.adjust.method = "fdr")
# ZT1    
# ZT13 6.8e-15

pairwise.t.test(rpca$PC1, rpca$ZT_nash, p.adjust.method = "fdr")

# ZT1_NASH ZT13_Non_NASH
# ZT13_Non_NASH 3.2e-13  -            
#   ZT13_NASH     3.7e-10  0.76   

plot_a<-plot_grid(NULL,cond_rpca, nrow=2,rel_heights = c(3, 10))
plot_b<-plot_grid(plot_a,p, rel_widths = c(1, 4))
#plot_c<-plot_grid(NULL,nash_rpca, rel_widths = c(1,5.5))
plot_c<-plot_grid(NULL,ZT_rpca, rel_widths = c(1,11))
final_plot <- plot_grid(plot_b,plot_c,nrow=2,rel_heights  = c(3, 1))

#ggsave("transcriptomics/rpca_results_liver/rpca_results_NASH_FAFTprot/SFR23_1212_NASH_liverRNA_NASHcat_wboxplot_RPCA_v2.pdf", plot=final_plot,height=5, width=5.7)
ggsave("transcriptomics/rpca_results_liver/rpca_results_NASH_FAFTprot/SFR24_0417_NASH_liverRNA_NASHZTcat_wboxplot_RPCA_v2.pdf", plot=final_plot,height=5, width=5.7)

##################################################################

#deseq2 results

#comparecluster with formula option

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
  #formula_res <- compareCluster(Entrez~direction+condition, OrgDb = "org.Mm.eg.db", data=mydf,
  #                              fun="enrichGO",pvalueCutoff=0.05,pAdjustMethod="fdr")
  formula_res <- compareCluster(Entrez~direction+condition, organism = "mmu", data=mydf,
                                fun="enrichKEGG",pvalueCutoff=0.05,pAdjustMethod="fdr")
  return(formula_res)
}

#nash_category

organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
require(DOSE)

deq_nash<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.NASH_category-NASH-Non_NASH-pval_05.csv")
gse_nash_plot(deq_nash)+ggtitle("Non-NASH vs NASH")+ scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=12),
        legend.position = "left")
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-NASH-Non_NASH.gse.pdf",height=6, width=10)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-NASH-Non_NASH.kegg.pdf",height=6, width=10)
# deq_nash_time<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.NASH_category-ZT1.NASH-ZT13.NASH-pval_05.csv")
# gse_nash_plot(deq_nash_time)+ggtitle("ZT1 NASH vs ZT13 NASH")
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-ZT1.NASH-ZT13.NASH.gse.pdf",height=8, width=7)

#deq_nash_time2<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.NASH_category-ZT13.NASH-ZT13.Non_NASH-pval_05.csv")
#gse_nash_plot(deq_nash_time2)+ggtitle("Non-NASH vs NASH ZT13") no hits

deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FA.Non_NASH.csv')%>%
  filter(pvalue<0.05)
deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.NASH-ZT13.FT.Non_NASH.csv')%>%
  filter(pvalue<0.05)

FAFT_nash<-gse_compareclust(deq_nash_FA,deq_nash_FT,"FA","FT")
dotplot(FAFT_nash,x="direction", showCategory=100) +
  facet_grid(~condition) + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(angle = 45, vjust = 0.5))
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-nashnonnash-FAFT.compclust_v2.pdf",height=10, width=7)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-nashnonnash-FAFT.compclust_v2.pdf",height=8, width=9)


# deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FA.Non_NASH-pval_05.csv')
deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FA.Non_NASH.csv')
gse_nash_plot(deq_nash_FA)+ggtitle("FA non-NASH vs NASH")+
  scale_y_discrete(position = "right") 
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FA.NASH-FA.Non_NASH.gse.pdf",height=6, width=12)
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FA.NASH-FA.Non_NASH.kegg.pdf",height=6, width=12)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT13FA.NASH-FA.Non_NASH.kegg.pdf",height=6, width=12)

# deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FT.NASH-FT.Non_NASH-pval_05.csv')
deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.NASH-ZT13.FT.Non_NASH.csv')
gse_nash_plot(deq_nash_FT)+ggtitle("FT non-NASH vs NASH")+
  scale_y_discrete(position = "right") 
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FT.NASH-FT.Non_NASH.gse.pdf",height=6, width=7)
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FT.NASH-FT.Non_NASH.kegg.pdf",height=6, width=12)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT13FT.NASH-FT.Non_NASH.kegg.pdf",height=6, width=12)


#plot together using comparecluster
#deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FA.Non_NASH-pval_05.csv')
#deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FT.NASH-FT.Non_NASH-pval_05.csv')
#just ZT13 as there's no nonmash for ZT1
deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FA.Non_NASH.csv')%>%
  filter(pvalue<0.05)
deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.NASH-ZT13.FT.Non_NASH.csv')%>%
  filter(pvalue<0.05)


FAFT_nash<-gse_compareclust(deq_nash_FA,deq_nash_FT,"FA","FT")
dotplot(FAFT_nash,x="direction", showCategory=100) +
  facet_grid(~condition) + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(angle = 45, vjust = 0.5))
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-nashnonnash-FAFT.compclust_v2.pdf",height=10, width=7)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-nashnonnash-FAFT.compclust_v2.pdf",height=8, width=9)

#plotting kegg
FAFT_nash<-gse_compareclust(deq_nash_FA,deq_nash_FT,"FA","FT")
FAFT_nash_df<-FAFT_nash@compareClusterResult#%>%
  #filter(!(category=="Human Diseases" | is.na(category)))
FAFT_nash@compareClusterResult = FAFT_nash@compareClusterResult[FAFT_nash@compareClusterResult$ID %in% FAFT_nash_df$ID, ]
dotplot(FAFT_nash,x="direction", showCategory=60) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(angle = 45, vjust = 0.5))
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-nashnonnash-FAFT.compclust_v2.kegg.pdf",height=10, width=10)
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-nashnonnash-FAFT.compclust_v2.keggsubset.pdf",height=10, width=10)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-FAFT.compclust_v2.kegg.pdf",height=9, width=10)

#plot together just the shared FA and FT hits using comparecluster
# deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FA.Non_NASH-pval_05.csv')%>%
#   filter(V1 %in% just_shared)
# deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FT.NASH-FT.Non_NASH-pval_05.csv')%>%
#   filter(V1 %in% just_shared)
deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FA.Non_NASH.csv')%>%
  filter(pvalue<0.05)%>%filter(V1 %in% just_shared)
deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.NASH-ZT13.FT.Non_NASH.csv')%>%
  filter(pvalue<0.05)%>%filter(V1 %in% just_shared) #note no functional enrichment found for overlap

FAFT_nash<-gse_compareclust(deq_nash_FA,deq_nash_FT,"FA","FT")
dotplot(FAFT_nash,x="direction", showCategory=15) +
  facet_grid(~condition) + scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left")
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-nashnonnash-FAFT.compclust_justshared.pdf",height=6, width=7)

FAFT_nash<-gse_compareclust(deq_nash_FA,deq_nash_FT,"FA","FT")
FAFT_nash_df<-FAFT_nash@compareClusterResult%>%
  filter(!(category=="Human Diseases" | is.na(category)))
FAFT_nash@compareClusterResult = FAFT_nash@compareClusterResult[FAFT_nash@compareClusterResult$ID %in% FAFT_nash_df$ID, ]

dotplot(FAFT_nash,x="direction", showCategory=15) +
  facet_grid(category~condition, scales="free",space='free_y') + scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left")
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-nashnonnash-FAFT.compclust_justshared.kegg.pdf",height=4, width=8)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-nashnonnash-FAFT.compclust_justshared_nohumandis.kegg.pdf",height=2.5, width=8)

#plot together just the non-overlap FA and FT hits using comparecluster
# deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FA.Non_NASH-pval_05.csv')%>%
#   filter(!(V1 %in% just_shared))
# deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FT.NASH-FT.Non_NASH-pval_05.csv')%>%
#   filter(!(V1 %in% just_shared))
deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FA.Non_NASH.csv')%>%
  filter(pvalue<0.05)%>%filter(!(V1 %in% just_shared))
deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.NASH-ZT13.FT.Non_NASH.csv')%>%
  filter(pvalue<0.05)%>%filter(!(V1 %in% just_shared))#note no functional enrichment found for overlap


FAFT_nash<-gse_compareclust(deq_nash_FA,deq_nash_FT,"FA","FT")
dotplot(FAFT_nash,x="direction", showCategory=15) +
  facet_grid(~condition) + scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left")

FAFT_nash<-gse_compareclust(deq_nash_FA,deq_nash_FT,"FA","FT")
FAFT_nash_df<-FAFT_nash@compareClusterResult#%>%
  # filter(!(category=="Human Diseases" | is.na(category)))
FAFT_nash@compareClusterResult = FAFT_nash@compareClusterResult[FAFT_nash@compareClusterResult$ID %in% FAFT_nash_df$ID, ]

dotplot(FAFT_nash,x="direction", showCategory=53) +
  facet_grid(category~condition, scales="free",space='free_y') + scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left")
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-nashnonnash-FAFT.compclust_justshared.kegg.pdf",height=4, width=8)
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-nashnonnash-FAFT.compclust_nonoverlap_nohumandis.kegg.pdf",height=10, width=10)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-FAFT.compclust_justshared.kegg.pdf",height=9, width=9)
#####

#get all the genes under non-nash from 193 shared
annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")

#get ppar genes
FAFT_nash_df<-FAFT_nash%>%as.data.frame()

nonNASH_ppar<-FAFT_nash_df%>%filter(direction=="a-Non-NASH" & grepl("PPAR signaling pathway",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="PPAR")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

#hits from kegg
nonNASH_carb<-FAFT_nash_df%>%filter(direction=="a-Non-NASH" & grepl("Glyoxylate and dicarboxylate",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="carb_mtb")

nonNASH_aa<-FAFT_nash_df%>%filter(direction=="a-Non-NASH" & grepl("Arginine and proline metabolism",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="aa_mtb")

NASH_immune<-FAFT_nash_df%>%filter(direction=="b-NASH" & grepl("RIG",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="immune")

sel_genes<-rbind(nonNASH_carb,nonNASH_aa,NASH_immune)%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)
md<-fread("transcriptomics/liver_metadata_wnashscore_fibste.tsv") %>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))

#plot heatmap of these specific hits (Z-score)
genecounts<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.vsd.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  #filter(gene_id %in% sel_genes$gene_id)%>%
  filter(gene_id %in% nonNASH_ppar$gene_id)%>%
  gather(sampleid,TPM_counts,-gene_id) %>%
  left_join(.,md,by="sampleid") %>%
  filter(condition!="NA")%>%
  left_join(.,nonNASH_ppar,by="gene_id")%>%
  group_by(gene_id,condition,NASH_category,Kegg)%>%dplyr::summarise(mn_TPM=mean(TPM_counts))%>%
  group_by(gene_id)%>%mutate(Zscore=(mn_TPM - mean(mn_TPM))/sd(mn_TPM))%>%
  mutate(condition=factor(condition,levels=c("FA","FT")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","not_applicable")))%>%
         #Kegg=factor(Kegg,levels=c("carb_mtb","aa_mtb","immune")))%>%
  separate(gene_id, c("gene_number","gene_symbol"), remove=FALSE)

# Run clustering
#genecounts_df<-fread("transcriptomics/data_files/liver_gene_tpm_matrix.csv")%>%
genecounts_df<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.vsd.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  #filter(gene_id %in% sel_genes$gene_id) 
  filter(gene_id %in% nonNASH_ppar$gene_id)
genecounts_matrix<-genecounts_df%>%column_to_rownames("gene_id") %>%as.matrix()
genecounts_dendro <- as.dendrogram(hclust(d = dist(x = genecounts_matrix)))
dendro_plot <- ggdendrogram(data = genecounts_dendro, rotate = TRUE)

genecounts_order <- order.dendrogram(genecounts_dendro)
genecounts$gene_id <- factor(x = genecounts$gene_id,
                             levels = genecounts_df$gene_id[genecounts_order], 
                             ordered = TRUE)


plt<-ggplot(genecounts,aes(x=NASH_category, y=gene_id)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+scale_y_discrete(position = "right") + 
  facet_grid(Kegg~condition,scales="free",space="free",switch = "y")+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),axis.text.y = element_text(size = 6),
        panel.spacing.y=unit(0.01, "lines"),strip.text.y = element_text(angle = 0))+
  scale_fill_viridis_b(option = "plasma")

#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_kegg_hits_GOenrich_shared.pdf", plot=plt,height=3, width=4.5)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_kegg_hits_GOenrich_PPAR.pdf", plot=plt,height=2.5, width=4.5)




#hits from 
nonNASH_lyase<-FAFT_nash_df%>%filter(direction=="a-Non-NASH" & grepl("lyase",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(GO="lyase")

nonNASH_bind<-FAFT_nash_df%>%filter(direction=="a-Non-NASH" & grepl("binding",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(GO="binding")%>%
  filter(gene_symbol!="Chac1")

NASH_ubi<-FAFT_nash_df%>%filter(direction=="b-NASH" & grepl("ubiquitin",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(GO="ubiquitin")


#get all the genes with ubiquitin or ubiquitin-like activity enriched under NASH FA and FT
annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")

nonNASH_ubi<-FAFT_nash_df%>%filter(Cluster=="b-NASH.FA" & grepl("ubiquitin",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(GO="ubiquitin")

#get all the genes with GTP activity enriched under NASH FA and FT

nonNASH_gtp<-FAFT_nash_df%>%filter(Cluster=="b-NASH.FT" & grepl("GTP",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(GO="GTP")

#get all the genes with ligase activity enriched under non-NASH FT

NASH_FT_ligase<-FAFT_nash_df%>%filter(Cluster=="a-Non-NASH.FT" & grepl("ligase",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(GO="ligase")

#get all the genes with lyase activity enriched under non-NASH FA
NASH_FA_lyase<-FAFT_nash_df%>%filter(Cluster=="a-Non-NASH.FA" & grepl("lyase",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(GO="lyase")

sel_genes<-rbind(nonNASH_lyase,nonNASH_bind,NASH_ubi)%>%
#sel_genes<-rbind(nonNASH_ubi,nonNASH_gtp)%>%
#sel_genes<-rbind(NASH_FT_ligase,NASH_FA_lyase)%>%
# sel_genes<-rbind(nonNASH_ubi)%>%
#   dplyr::rename(gene_id=stringtie_id)
#sel_genes<-rbind(nonNASH_gtp)%>%
  dplyr::rename(gene_id=stringtie_id)
md<-fread("transcriptomics/liver_metadata_wnashscore_fibste.tsv") %>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))

#plot heatmap of these specific hits (log2foldchange)

##ligase/lyase
deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FA.Non_NASH.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")
deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FT.NASH-FT.Non_NASH.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")

lyl_FA<-deq_nash_FA%>%filter(gene_symbol %in% sel_genes$gene_symbol)%>%mutate(condition="FA")%>%
  mutate(GO=ifelse(gene_symbol %in% NASH_FA_lyase$gene_symbol, "lyase","ligase" ))
lyl_FT<-deq_nash_FT%>%filter(gene_symbol %in% sel_genes$gene_symbol)%>%mutate(condition="FT")%>%
  mutate(GO=ifelse(gene_symbol %in% NASH_FT_ligase$gene_symbol, "ligase","lyase" ))

lyl<-rbind(lyl_FA,lyl_FT)%>%
  mutate(GO=factor(GO,levels=c("lyase","ligase")))

plt<-ggplot(lyl,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(aes(fill=log2FoldChange))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(rows = vars(GO),scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),axis.text.y = element_text(size = 7),
        panel.spacing.y=unit(0.01, "lines"),strip.text.y = element_text(angle = 0))+
  scale_fill_viridis_c(option = "plasma")

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-ligaselyase.pdf", plot=plt,height=4, width=3.5)

##ubiquitin/GTP
ubi_FA<-deq_nash_FA%>%filter(gene_symbol %in% sel_genes$gene_symbol)%>%mutate(condition="FA")%>%
  mutate(GO=ifelse(gene_symbol %in% nonNASH_ubi$gene_symbol, "ubiquitin","GTP" ))
gtp_FT<-deq_nash_FT%>%filter(gene_symbol %in% sel_genes$gene_symbol)%>%mutate(condition="FT")%>%
  mutate(GO=ifelse(gene_symbol %in% nonNASH_gtp$gene_symbol, "GTP","ubiquitin" ))

##ubiquitin
ubi_FA<-deq_nash_FA%>%filter(gene_symbol %in% nonNASH_ubi$gene_symbol)%>%mutate(condition="FA")
ubi_FT<-deq_nash_FT%>%filter(gene_symbol %in% nonNASH_ubi$gene_symbol)%>%mutate(condition="FT")

ubi<-rbind(ubi_FA,ubi_FT)

plt<-ggplot(ubi,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(aes(fill=log2FoldChange))+
  scale_x_discrete(expand = c(0, 0))+
  theme(axis.ticks.y=element_blank(),axis.text.y = element_text(size = 7),
        panel.spacing.y=unit(0.01, "lines"),strip.text.y = element_text(angle = 0))+
  scale_fill_viridis_c(option = "plasma")+ggtitle("ubiquitin &\nubiquitin-like activity")

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-ubiquitin.pdf", plot=plt,height=6, width=3.5)

##GTP
gtp_FA<-deq_nash_FA%>%filter(gene_symbol %in% nonNASH_gtp$gene_symbol)%>%mutate(condition="FA")
gtp_FT<-deq_nash_FT%>%filter(gene_symbol %in% nonNASH_gtp$gene_symbol)%>%mutate(condition="FT")

gtp<-rbind(gtp_FA,gtp_FT)

plt<-ggplot(gtp,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(aes(fill=log2FoldChange))+
  scale_x_discrete(expand = c(0, 0))+
  theme(axis.ticks.y=element_blank(),axis.text.y = element_text(size = 7),
        panel.spacing.y=unit(0.01, "lines"),strip.text.y = element_text(angle = 0))+
  scale_fill_viridis_c(option = "plasma")+ggtitle("GTP activity")

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-gtp.pdf", plot=plt,height=8, width=3.5)


#plot heatmap of these specific hits (Z-score)
#genecounts<-fread("transcriptomics/data_files/liver_gene_tpm_matrix.csv")
genecounts<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.vsd.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  filter(gene_id %in% sel_genes$gene_id)%>%
  gather(sampleid,TPM_counts,-gene_id) %>%
  left_join(.,md,by="sampleid") %>%
  filter(condition!="NA")%>%
  left_join(.,sel_genes,by="gene_id")%>%
  group_by(gene_id,condition,NASH_category,GO)%>%dplyr::summarise(mn_TPM=mean(TPM_counts))%>%
  group_by(gene_id)%>%mutate(Zscore=(mn_TPM - mean(mn_TPM))/sd(mn_TPM))%>%
  mutate(condition=factor(condition,levels=c("FA","FT")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","not_applicable")),
         #GO=factor(GO,levels=c("lyase","ligase","ubiquitin","GTP")))%>% #"ubiquitin",
         GO=factor(GO,levels=c("lyase","binding","ubiquitin")))%>% #"ubiquitin",
  separate(gene_id, c("gene_number","gene_symbol"), remove=FALSE)

# Run clustering
#genecounts_df<-fread("transcriptomics/data_files/liver_gene_tpm_matrix.csv")%>%
genecounts_df<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.vsd.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  filter(gene_id %in% sel_genes$gene_id) 
genecounts_matrix<-genecounts_df%>%column_to_rownames("gene_id") %>%as.matrix()
genecounts_dendro <- as.dendrogram(hclust(d = dist(x = genecounts_matrix)))
dendro_plot <- ggdendrogram(data = genecounts_dendro, rotate = TRUE)

genecounts_order <- order.dendrogram(genecounts_dendro)
genecounts$gene_id <- factor(x = genecounts$gene_id,
                               levels = genecounts_df$gene_id[genecounts_order], 
                               ordered = TRUE)


plt<-ggplot(genecounts,aes(x=NASH_category, y=gene_symbol)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+scale_y_discrete(position = "right") + 
  facet_grid(GO~condition,scales="free",space="free",switch = "y")+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),axis.text.y = element_text(size = 6),
        panel.spacing.y=unit(0.01, "lines"),strip.text.y = element_text(angle = 0))+
  scale_fill_viridis_b(option = "plasma")

#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_ligaselyase_hits_GOenrich_v5.pdf", plot=plt,height=4, width=4.5)
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_ubiquitin_hits_GOenrich.pdf", plot=plt,height=6, width=5)
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_gtp_hits_GOenrich.pdf", plot=plt,height=8, width=5)
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_ubigtp_hits_GOenrich_v2.pdf", plot=plt,height=10, width=4)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_lyasebindubi_hits_GOenrich_shared.pdf", plot=plt,height=4, width=4)


cnetplot(FAFT_nash)+scale_fill_manual(values=c("#E69F00","limegreen","#D55E00","darkgreen"))
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-nashnonnash-FAFT-network.compclust.pdf",height=15, width=15)

#log(FA.NASH/FT.NASH)
deq_nash_nash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FT.NASH-pval_05.csv')
gse_nash_plot(deq_nash_nash)+ggtitle("FT vs. FA NASH")+scale_y_discrete(position = "right") 
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FA.NASH-FT.NASH.gse.pdf",height=6, width=8)

#log(FA.non-NASH/FT.non-NASH)
deq_nash_nonnash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.Non_NASH-FT.Non_NASH-pval_05.csv')
gse_nash_plot(deq_nash_nonnash)+ggtitle("FT vs. FA non-NASH")+scale_y_discrete(position = "right") 
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FA.Non_NASH-FT.Non_NASH.gse.pdf",height=6, width=7)

deq_nash_nash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FT.NASH-pval_05.csv')
deq_nash_nonnash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.Non_NASH-FT.Non_NASH-pval_05.csv')
nonnash_nash<-gse_compareclust(deq_nash_nonnash,deq_nash_nash,"a-Non-NASH","b-NASH")

dotplot(nonnash_nash,x="direction", showCategory=15) + facet_grid(~condition)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FAFT-nashnonnash.compclust.pdf",height=10, width=10)

cnetplot(nonnash_nash)+scale_fill_manual(values=c("#E69F00","#D55E00","limegreen","darkgreen"))
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FAFT-nashnonnash-network.compclust.pdf",height=15, width=15)

#deq_steatosis_grd<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.steatosis_grade_new-0-1+-pval_05.csv') not enough samples
#deq_fibrosis_stg<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.fibrosis_stage_new-0-1+-pval_05.csv") not enough samples

# deq_timept<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-NA.ZT1-NA.ZT13-pval_05.csv")
# gse_nash_plot(deq_timept)+ggtitle("ZT1 vs ZT13 NA")
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint-NA.ZT1-NA.ZT13.gse.pdf",height=8, width=7)
# 
# deq_timept_FA<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-FA.ZT1-FA.ZT13-pval_05.csv")
# gse_nash_plot(deq_timept_FA)+ggtitle("ZT1 vs ZT13 FA")
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint-FA.ZT1-FA.ZT13.gse.pdf",height=8, width=7)
# 
# deq_timept_FT<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-FT.ZT1-FT.ZT13-pval_05.csv")
# gse_nash_plot(deq_timept_FT)+ggtitle("ZT1 vs ZT13 FT")
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint-FT.ZT1-FT.ZT13.gse.pdf",height=8, width=7)
# 
# deq_timept_FA<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-FA.ZT1-FA.ZT13-pval_05.csv")
# deq_timept_FT<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-FT.ZT1-FT.ZT13-pval_05.csv")
# ZT_FAFT<-gse_compareclust(deq_timept_FA,deq_timept_FA,"FA","FT")
# 
# dotplot(ZT_FAFT,x="direction", showCategory=25) + facet_grid(~condition)
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FT.ZT1-FT.ZT13.compclust.pdf",height=10, width=10)
# 
# cnetplot(ZT_FAFT)+scale_fill_manual(values=c("#E69F00","limegreen","#D55E00","darkgreen"))
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FT.ZT1-FT.ZT13-network.compclust.pdf",height=15, width=15)
##################################################################
#venn diagrams
# deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FA.Non_NASH-pval_05.csv')#%>%
#   #filter(log2FoldChange>0)
# deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FT.NASH-FT.Non_NASH-pval_05.csv')#%>%
#   #filter(log2FoldChange>0)
deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FA.Non_NASH.csv')%>%
  filter(pvalue<0.05)
  #filter(padj<0.05)
deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.NASH-ZT13.FT.Non_NASH.csv')%>%
  filter(pvalue<0.05)
  #filter(padj<0.05)
list_venn <- list(FA = deq_nash_FA$V1,
                  FT = deq_nash_FT$V1)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

just_shared<-all$`FA:FT`

p<-venn.diagram(list_venn, fill = c("#D55E00","#009E73"),height = 10,lty = 0, 
                main="ZT13 nonnash vs nash",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
#pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_venn_nonnash_nash.pdf")
#pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_venn_non-nash.pdf")
#pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_venn_nash.pdf")
pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_ZT13venn_nash.pdf")
#pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_ZT13venn_nash_padj.pdf")
grid.draw(p)
dev.off()

deq_nash_nash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FT.NASH-pval_05.csv')
deq_nash_nonnash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.Non_NASH-FT.Non_NASH-pval_05.csv')

list_venn <- list(NASH = deq_nash_nash$V1,
                  NonNASH = deq_nash_nonnash$V1)

p<-venn.diagram(list_venn, fill = c("#D55E00","#009E73"),height = 10,lty = 0, 
                main="FA vs. FT",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_venn_FAFT_nonnash_nash.pdf")
grid.draw(p)
dev.off()

# deq_timept_FA<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-FA.ZT1-FA.ZT13-pval_05.csv")
# deq_timept_FT<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-FT.ZT1-FT.ZT13-pval_05.csv")
# 
# list_venn <- list(FA = deq_timept_FA$V1,
#                   FT = deq_timept_FT$V1)
# 
# p<-venn.diagram(list_venn, fill = c("#D55E00","#009E73"),height = 10,lty = 0, 
#                 main="ZT1 vs ZT13",
#                 width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)
# 
# grid.draw(p)
# pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_venn_ZT1ZT13.pdf")
# grid.draw(p)
# dev.off()
##################################################################

#make volcano plots
annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")%>%
  filter(molecule_type=="protein_coding")

#deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FA.Non_NASH.csv')%>%
deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FA.Non_NASH.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%
  left_join(.,annot,by="stringtie_id")

# p<-EnhancedVolcano(
#   deq_nash_FA,
#   lab = deq_nash_FA$gene_symbol,
#   title = 'non-NASH vs. NASH (FA)',
#   x = "log2FoldChange",
#   y = "padj",
#   labSize = 2.0,
#   pointSize = 1,
#   col=c('gray37', 'gray37', 'darkorange3', 'darkgoldenrod2'),
#   pCutoff = 0.05
# ) +theme_pubr()+theme(text=element_text(size=7))
# 
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FA.NASH-FA.Non_NASH.volplot.pdf",height=3, width=3)

p<-EnhancedVolcano(
  deq_nash_FA,
  lab = deq_nash_FA$gene_symbol,
  xlim=c(-5,5),
  ylim=c(0,6),
  title = 'non-NASH vs. NASH (FA)',
  subtitle=NA,
  caption=NA,
  x = "log2FoldChange",
  y = "pvalue",
  labSize = 0,
  pointSize = 1,
  colAlpha=1,
  col=c('gray37', 'gray37','darkgoldenrod2','darkorange3'),
  pCutoff = 0.05,
  legendPosition = "none",
  border="full",
  gridlines.major=FALSE,
  gridlines.minor=FALSE
) +theme(text=element_text(size=7))

# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FA.NASH-FA.Non_NASH.volplot_v2.pdf",height=5, width=5)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-FA.NASH-FA.Non_NASH.volplot_v2.pdf",height=5, width=5)


#deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FT.NASH-FT.Non_NASH.csv')%>%
deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.NASH-ZT13.FT.Non_NASH.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%
  left_join(.,annot,by="stringtie_id")

# p<-EnhancedVolcano(
#   deq_nash_FT,
#   lab = deq_nash_FT$gene_symbol,
#   title = 'non-NASH vs. NASH (FT)',
#   x = "log2FoldChange",
#   y = "padj",
#   labSize = 2.0,
#   pointSize = 1,
#   col=c('gray37', 'gray37', "darkgreen","chartreuse3"),
#   pCutoff = 0.05
# ) +theme_pubr()+theme(text=element_text(size=7))
# 
# ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FT.NASH-FT.Non_NASH.volplot.pdf",height=3, width=3)

p<-EnhancedVolcano(
  deq_nash_FT,
  lab = deq_nash_FT$gene_symbol,
  xlim=c(-5,5),
  ylim=c(0,6),
  title = 'non-NASH vs. NASH (FT)',
  subtitle=NA,
  caption=NA,
  x = "log2FoldChange",
  y = "pvalue",
  labSize = 0,
  colAlpha=1,
  pointSize = 1,
  col=c('gray37', 'gray37',"chartreuse3","darkgreen"),
  pCutoff = 0.05,
  legendPosition = "none",
  border="full",
  gridlines.major=FALSE,
  gridlines.minor=FALSE
) +theme(text=element_text(size=7))

#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-FT.NASH-FT.Non_NASH.volplot_v2.pdf",height=5, width=5)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-FT.NASH-FT.Non_NASH.volplot_v2.pdf",height=5, width=5)


deq_nash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.NASH_category-NASH-Non_NASH.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%
  left_join(.,annot,by="stringtie_id")

p<-EnhancedVolcano(
  deq_nash,
  lab = deq_nash$gene_symbol,
  title = 'non-NASH vs. NASH',
  x = "log2FoldChange",
  y = "padj",
  labSize = 2.0,
  pointSize = 1,
  pCutoff = 0.05
) +theme_pubr()+theme(text=element_text(size=7))

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.NASH_category-NASH-Non_NASH.volplot.pdf",height=3, width=3)

##################################################################

#make heatmap of functions of interest
clock_genes<-c("Npas2","Arntl","Clock","Nr0b2","Nr1d1","Nr1d2","Per1","Rora")
BA_genes<-c("Cyp27a1","Cyp7a1","Cyp7b1","Cyp8b1","Nr1h4","Osbpl9")
annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")

deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FA.Non_NASH.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")
deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FT.NASH-FT.Non_NASH.csv')%>%
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
                       limits=c(-7,7))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  geom_text(aes(label = signif(padj,digits=2)), color = "black") +
  guides(fill = guide_colourbar(barheight = 0.5))+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0), legend.positio="top") + 
 labs(x="",y="")+ggtitle("Circadian Clock Genes\nNon-NASH v. NASH") # scale_fill_viridis_c() 

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-clockgenes_v2.pdf", plot=plt,height=4, width=3)


BA_FA<-deq_nash_FA%>%filter(gene_symbol %in% BA_genes)%>%mutate(condition="FA")
BA_FT<-deq_nash_FT%>%filter(gene_symbol %in% BA_genes)%>%mutate(condition="FT")

BA<-rbind(BA_FA,BA_FT)

plt<-ggplot(BA,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(color="black", aes(fill=log2FoldChange))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000", #colors in the scale
                       midpoint=0,    #same midpoint for plots (mean of the range)
                       breaks=seq(-100,100,0.5), #breaks in the scale bar
                       limits=c(-1.5,1.5))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  geom_text(aes(label = signif(padj,digits=2)), color = "black") +
  guides(fill = guide_colourbar(barheight = 0.5))+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0),legend.position = "top") + 
  labs(x="",y="")+ggtitle("Bile Acid Genes\nNon-NASH v. NASH")

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-BAgenes_v2.pdf", plot=plt,height=4, width=3)


deq_timept_FA<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-FA.ZT1-FA.ZT13.csv")%>%
  dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")
deq_timept_FT<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-FT.ZT1-FT.ZT13.csv")%>%
  dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")
deq_timept_NA<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-NA.ZT1-NA.ZT13.csv")%>%
  dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")

clk_FA<-deq_timept_FA%>%filter(gene_symbol %in% clock_genes)%>%mutate(condition="FA")
clk_FT<-deq_timept_FT%>%filter(gene_symbol %in% clock_genes)%>%mutate(condition="FT")
clk_NA<-deq_timept_NA%>%filter(gene_symbol %in% clock_genes)%>%mutate(condition="NA")

clk<-rbind(clk_FA,clk_FT,clk_NA)%>%
  mutate(log2FoldChange=-1*log2FoldChange,
         condition=factor(condition,levels=c("NA","FA","FT")))

plt<-ggplot(clk,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(color="black", aes(fill=log2FoldChange))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000", #colors in the scale
                       midpoint=0,    #same midpoint for plots (mean of the range)
                       breaks=seq(-100,100,2), #breaks in the scale bar
                       limits=c(-7,7))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  geom_text(aes(label = signif(padj,digits=2)), color = "black") +
  guides(fill = guide_colourbar(barheight = 0.5))+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0),legend.position = "top") +
  labs(x="",y="")+ggtitle("Circadian Clock Genes\nZT1 vs. ZT13")

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-clockgenes_v2.pdf", plot=plt,height=4, width=3.5)

BA_FA<-deq_timept_FA%>%filter(gene_symbol %in% BA_genes)%>%mutate(condition="FA")
BA_FT<-deq_timept_FT%>%filter(gene_symbol %in% BA_genes)%>%mutate(condition="FT")
BA_NA<-deq_timept_NA%>%filter(gene_symbol %in% BA_genes)%>%mutate(condition="NA")

BA<-rbind(BA_FA,BA_FT,BA_NA)%>%
  mutate(log2FoldChange=-1*log2FoldChange,
         condition=factor(condition,levels=c("NA","FA","FT")))

plt<-ggplot(BA,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(color="black", aes(fill=log2FoldChange))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000", #colors in the scale
                       midpoint=0,    #same midpoint for plots (mean of the range)
                       breaks=seq(-100,100,0.5), #breaks in the scale bar
                       limits=c(-1.5,1.5))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  geom_text(aes(label = signif(padj,digits=2)), color = "black") +
  guides(fill = guide_colourbar(barheight = 0.5))+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0),legend.position = "top") +
  labs(x="",y="")+ggtitle("Bile Acid Genes\nZT1 vs. ZT13")

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-BAgenes_v2.pdf", plot=plt,height=4, width=3.5)


