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
library(UpSetR)
###############################################################
#NASH all (prot only)
ord <- read_qza("transcriptomics/rpca_results_liver/rpca_results_NASH_prot/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"transcriptomics/rpca_results_liver/rpca_results_NASH_prot/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"transcriptomics/rpca_results_liver/rpca_results_NASH_prot/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("transcriptomics/liver_metadata_wnashscore_fibste.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         NASH_category= ifelse(NASH_category=="not_applicable", "Non_NASH",NASH_category))%>%
  dplyr::rename(sample_name=sampleid)%>%
  mutate(mashZT=paste(timepoint,NASH_category,sep="_"))
write.table(md,"transcriptomics/liver_metadata_wnashscore_fibste_wmashZT.tsv",sep = "\t",row.names = FALSE,quote=FALSE)


rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(phase=ifelse(timepoint=="ZT1","light","dark"))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH")),
         mashZT=factor(mashZT,levels=c("ZT1_Non_NASH","ZT1_NASH","ZT13_Non_NASH","ZT13_NASH")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, shape=mashZT)) +
  geom_point(aes(fill=condition),alpha=1.0, size=2.5) +
  theme_pubr() +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_shape_manual(values=c(21,22,23,24)) +
  labs(x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("NASH liver RNA (Wk 12)")+ theme(plot.title = element_text(face = "bold"))
ggsave("transcriptomics/rpca_results_liver/rpca_results_NASH_prot/SFR24_0419_NASH_liverRNA_NASHcat_RPCA.pdf", plot=p,height=4, width=4)

##add box plot to axes
cond_rpca <- ggplot(rpca, aes(x =factor(condition, levels = rev(levels(condition))), y = PC1,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_point(colour="black", aes(fill=condition),
                                       position=position_jitter(width=0.1,height = 0.1),
                                       size=3,pch=21) + 
  labs(x = NA, y = "PC1") +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  # theme_void()
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")+coord_flip()

pairwise.t.test(rpca$PC1, rpca$condition, p.adjust.method = "fdr")

# NA      FA     
# FA 8.3e-14 -      
#   FT < 2e-16 0.00016

ZT_rpca <- ggplot(rpca, aes(x =mashZT, y = PC2,shape=mashZT)) +
  geom_boxplot(alpha=0.3) + geom_point(fill="black",colour="white",position=position_jitter(width=0.1,height = 0.1),size=3) +
  scale_shape_manual(values=c(21,22,23,24)) +
  labs(x = NA, y = "PC2") +
  # theme_void()
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

pairwise.t.test(rpca$PC2, rpca$mashZT, p.adjust.method = "fdr")

# ZT1_Non_NASH ZT1_NASH ZT13_Non_NASH
# ZT1_NASH      0.0054       -        -            
#   ZT13_Non_NASH 2.3e-16      2.3e-16  -            
#   ZT13_NASH     9.3e-14      8.6e-13  0.0557       

plot_a<-plot_grid(NULL,ZT_rpca, nrow=2,rel_heights = c(2, 9))
plot_b<-plot_grid(plot_a,p, rel_widths = c(1, 3))
plot_c<-plot_grid(NULL,cond_rpca, rel_widths = c(2,5))
final_plot <- plot_grid(plot_b,plot_c,nrow=2,rel_heights  = c(3, 1))

ggsave("transcriptomics/rpca_results_liver/rpca_results_NASH_prot/SFR24_0419_NASH_liverRNA_NASHcat_wboxplot_RPCA.pdf", plot=final_plot,height=5, width=5.7)


###############################################################

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
  kk2 <- gseKEGG(geneList     = gene_list,
                 organism     = "mmu",
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 keyType       = "ncbi-geneid")

  p<-dotplot(kk2, showCategory=15, split=".sign") + facet_grid(.~.sign)
  return(p)
}

gse_compareclust<-function(df1,df2,df3,cond1,cond2,cond3){
#gse_compareclust<-function(df1,cond1){ 
#gse_compareclust<-function(df1,df2,cond1,cond2){
  gene_list1<-get_genelist(df1)
  gene_list2<-get_genelist(df2)
  gene_list3<-get_genelist(df3)
  
  mydf1 <- data.frame(Entrez=names(gene_list1), FC=gene_list1)%>%
    mutate(direction=ifelse(FC<0,"FA","FT"),
           condition=cond1)
  
  mydf2 <- data.frame(Entrez=names(gene_list2), FC=gene_list2)%>%
     mutate(direction=ifelse(FC<0,"FA","FT"),
            condition=cond2)
  
  mydf3 <- data.frame(Entrez=names(gene_list3), FC=gene_list3)%>%
     mutate(direction=ifelse(FC<0,"FA","FT"),
            condition=cond3)

  mydf<-rbind(mydf1,mydf2,mydf3)
  #mydf<-rbind(mydf1,mydf2)
  #formula_res <- compareCluster(Entrez~direction+condition, OrgDb = "org.Mm.eg.db", data=mydf,
  #                              fun="enrichGO",pvalueCutoff=0.05,pAdjustMethod="fdr")
  formula_res <- compareCluster(Entrez~direction+condition, organism = "mmu", data=mydf,
                                fun="enrichKEGG",pvalueCutoff=0.05,pAdjustMethod="fdr")
  return(formula_res)
}


gse_compareclust<-function(df1,df2,df3,cond1,cond2,cond3){
  #gse_compareclust<-function(df1,cond1){ 
  #gse_compareclust<-function(df1,df2,cond1,cond2){
  gene_list1<-get_genelist(df1)
  gene_list2<-get_genelist(df2)
  gene_list3<-get_genelist(df3)
  
  mydf1 <- data.frame(Entrez=names(gene_list1), FC=gene_list1)%>%
    mutate(direction=ifelse(FC<0,"FA","FT"),
           condition=cond1)
  
  mydf2 <- data.frame(Entrez=names(gene_list2), FC=gene_list2)%>%
    mutate(direction=ifelse(FC<0,"FA","FT"),
           condition=cond2)
  
  mydf3 <- data.frame(Entrez=names(gene_list3), FC=gene_list3)%>%
    mutate(direction=ifelse(FC<0,"FA","FT"),
           condition=cond3)
  
  mydf<-rbind(mydf1,mydf2,mydf3)
  #mydf<-rbind(mydf1,mydf2)
  #formula_res <- compareCluster(Entrez~direction+condition, OrgDb = "org.Mm.eg.db", data=mydf,
  #                              fun="enrichGO",pvalueCutoff=0.05,pAdjustMethod="fdr")
  formula_res <- compareCluster(Entrez~direction+condition, organism = "mmu", data=mydf,
                                fun="enrichKEGG",pvalueCutoff=0.05,pAdjustMethod="fdr")
  return(formula_res)
}

gse_compareclust_ZT<-function(df1,cond1){
  gene_list1<-get_genelist(df1)

  mydf1 <- data.frame(Entrez=names(gene_list1), FC=gene_list1)%>%
    mutate(direction=ifelse(FC<0,"a-ZT13_non-MASH","b-ZT1-MASH"),
           condition=cond1)

  formula_res <- compareCluster(Entrez~direction+condition, organism = "mmu", data=mydf1,
                                fun="enrichKEGG",pvalueCutoff=0.05,pAdjustMethod="fdr")
  return(formula_res)
}

gse_compareclust_allcond<-function(df1,df2,df3,cond1,cond2,cond3){
  #gse_compareclust_allcond<-function(df1,df3,cond1,cond3){
  gene_list1<-get_genelist(df1)
  gene_list2<-get_genelist(df2)
  gene_list3<-get_genelist(df3)
  
  mydf1 <- data.frame(Entrez=names(gene_list1), FC=gene_list1)%>%
    mutate(direction=ifelse(FC<0,"a-NA","b-FA"),
           condition=cond1)
  
  mydf2 <- data.frame(Entrez=names(gene_list2), FC=gene_list2)%>%
    mutate(direction=ifelse(FC<0,"a-NA","c-FT"),
           condition=cond2)

  mydf3 <- data.frame(Entrez=names(gene_list3), FC=gene_list3)%>%
    mutate(direction=ifelse(FC<0,"c-FT","b-FA"),
           condition=cond3)
  mydf<-rbind(mydf1,mydf2,mydf3)
  #mydf<-rbind(mydf1,mydf3)
  formula_res <- compareCluster(Entrez~direction+condition, organism = "mmu", data=mydf,
                                fun="enrichKEGG",pvalueCutoff=0.05,pAdjustMethod="fdr")
  return(formula_res)
}

###############################################################

deq_nash_ZT113FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-FA.ZT1-FA.ZT13-pval_05.csv')
deq_nash_ZT113FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-FT.ZT1-FT.ZT13-pval_05.csv')

deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FA.Non_NASH-pval_05.csv')
deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FT.NASH-FT.Non_NASH-pval_05.csv')

list_venn <- list(ZT1vZT13FA = deq_nash_ZT113FA$V1,
                  nonmashmashFA = deq_nash_FA$V1)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

justnonmashmash_FA<-all$nonmashmashFA

p<-venn.diagram(list_venn, fill = c('darkgoldenrod2','darkorange3'),height = 10,lty = 0, 
                main="FA time vs. disease",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_ZT1ZT13_nonmashmash_FA.pdf")
grid.draw(p)
dev.off()


list_venn <- list(ZT1vZT13FT = deq_nash_ZT113FT$V1,
                  nonmashmashFT = deq_nash_FT$V1)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

justnonmashmash_FT<-all$nonmashmashFT

p<-venn.diagram(list_venn, fill = c("chartreuse3","darkgreen"),height = 10,lty = 0, 
                main="FT time vs. disease",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_ZT1ZT13_nonmashmash_FT.pdf")
grid.draw(p)
dev.off()

deq_nash_FA<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FA.NASH-FA.Non_NASH-pval_05.csv')%>%
  filter(V1 %in% justnonmashmash_FA)
deq_nash_FT<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-FT.NASH-FT.Non_NASH-pval_05.csv')%>%
  filter(V1 %in% justnonmashmash_FT)

list_venn <- list(FA = deq_nash_FA$V1,
                  FT = deq_nash_FT$V1)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

p<-venn.diagram(list_venn, fill = c("#D55E00","#009E73"),height = 10,lty = 0, 
                main="FT time vs. disease",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_nonmashmashFAFT_noZT.pdf")
grid.draw(p)
dev.off()

FAFT_nash<-gse_compareclust(deq_nash_FA,deq_nash_FT,"FA","FT")
dotplot(FAFT_nash,x="direction", showCategory=200) +
  facet_grid(~condition) + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5)) #not much to show

justFA<-deq_nash_FA%>%filter(V1 %in% justnonmashmash_FA)%>%mutate(condition="FA")
justFT<-deq_nash_FT%>%filter(V1 %in% justnonmashmash_FT)%>%mutate(condition="FT")

justFAFT<-rbind(justFA,justFT)

md<-fread("transcriptomics/liver_metadata_wnashscore_fibste.tsv") %>%
  mutate(condition=ifelse(is.na(condition),"NA", condition))%>%
  mutate(mashZT=paste(timepoint,NASH_category,sep="_"))

#plot heatmap of these specific hits (Z-score)
#genecounts<-fread("transcriptomics/data_files/liver_gene_tpm_matrix.csv")
genecounts<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.vsd.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  filter(gene_id %in% justnonmashmash_FA)%>%
  gather(sampleid,TPM_counts,-gene_id) %>%
  left_join(.,md,by="sampleid") %>%
  filter(condition!="NA")%>%
  filter(mashZT!="ZT13_NASH")%>%
  #left_join(.,sel_genes,by="gene_id")%>%
  group_by(gene_id,condition,mashZT)%>%dplyr::summarise(mn_TPM=mean(TPM_counts))%>%
  group_by(gene_id)%>%mutate(Zscore=(mn_TPM - mean(mn_TPM))/sd(mn_TPM))%>%
  mutate(condition=factor(condition,levels=c("FA","FT")),
         mashZT=factor(mashZT,levels=c("ZT1_NASH","ZT13_Non_NASH")))%>%
         #GO=factor(GO,levels=c("lyase","ligase","ubiquitin","GTP")))%>% #"ubiquitin",
         #GO=factor(GO,levels=c("lyase","binding","ubiquitin")))%>% #"ubiquitin",
  separate(gene_id, c("gene_number","gene_symbol"), remove=FALSE)

# Run clustering
#genecounts_df<-fread("transcriptomics/data_files/liver_gene_tpm_matrix.csv")%>%
genecounts_df<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.vsd.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  filter(gene_id %in% justnonmashmash_FA) 
genecounts_matrix<-genecounts_df%>%column_to_rownames("gene_id") %>%as.matrix()
genecounts_dendro <- as.dendrogram(hclust(d = dist(x = genecounts_matrix)))
dendro_plot <- ggdendrogram(data = genecounts_dendro, rotate = TRUE)

genecounts_order <- order.dendrogram(genecounts_dendro)
genecounts$gene_id <- factor(x = genecounts$gene_id,
                             levels = genecounts_df$gene_id[genecounts_order], 
                             ordered = TRUE)

plt<-ggplot(genecounts,aes(x=mashZT, y=gene_symbol)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+scale_y_discrete(position = "right") + 
  facet_grid(~condition,scales="free",space="free",switch = "y")+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),axis.text.y = element_text(size = 6),
        panel.spacing.y=unit(0.01, "lines"),strip.text.y = element_text(angle = 0))+
  scale_fill_viridis_b(option = "plasma")

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.NASH_category-justFAnoZT.pdf", plot=plt,height=8, width=4)


###############################################################
#plot just the insulin related results from ZT1mash vs ZT13 nonmash
ZT1mash_ZT13nonmash<-fread('transcriptomics/NASH_liver_deseq2/other_comparisons/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT1.FT.NASH-ZT13.FT.Non_NASH.csv')%>%
  filter(V1 %in% nonNASH_insulin$gene_id)%>%
  dplyr::rename(gene_id=V1)%>%
  separate(gene_id, c("gene_number","gene_symbol"), remove=FALSE)%>%
  mutate(significance=ifelse(padj<0.05, "sig","notsig"),
         gene_symbol = fct_reorder(gene_symbol,log2FoldChange))

ggplot(ZT1mash_ZT13nonmash, aes(x =gene_symbol , y = log2FoldChange,fill=significance)) + 
  geom_bar(stat="identity", alpha=0.7)+ theme_pubr()+
  labs(title="Insulin Signaling Pathway",y="log2(FT_ZT1_MASH/FT_ZT13_non-MASH)")+
  scale_fill_manual(values=c("gray60","red4"))+
  scale_y_continuous(limits=c(-2,2))+
  theme(legend.position = "top")+
  coord_flip()
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1.FT.NASH-ZT13.FT.Non_NASH.insulin.pdf",height=6, width=4)

###############################################################
ZT1mash_ZT13nonmash<-fread('transcriptomics/NASH_liver_deseq2/other_comparisons/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT1.FT.NASH-ZT13.FT.Non_NASH-pval_05.csv')%>%
  #filter(V1 %in% nonNASHFT_insulin$gene_id)
  filter(V1 %in% nonNASHFT_ampk$gene_id)
gse_nash_plot(deq_nash_FANA_ZT13NAnomash)+ggtitle("FT ZT1 MASH vs. ZT13 non-MASH")+
  scale_y_discrete(position = "right") 
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1.FT.NASH-ZT13.FT.Non_NASH.kegg.pdf",height=6, width=12)

ZT_nash<-gse_compareclust_ZT(ZT1mash_ZT13nonmash,"ZT13 non-MASH vs. ZT1 MASH")
ZT_nash_df<-ZT_nash@compareClusterResult%>%
  filter(!(category=="Human Diseases" | is.na(category)))
ZT_nash@compareClusterResult = ZT_nash@compareClusterResult[ZT_nash@compareClusterResult$ID %in% ZT_nash_df$ID, ]
dotplot(ZT_nash,x="direction", showCategory=30) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1.FT.NASH-ZT13.FT.Non_NASH.nohumandisease.kegg.pdf",height=6, width=7.5)



deq_nash_FANA_ZT13NAnomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.NA.not_applicable-pval_05.csv')
deq_nash_FTNA_ZT13NAnomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.Non_NASH-ZT13.NA.not_applicable-pval_05.csv')
deq_nash_ZT113_NAnonmash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT1.NA.not_applicable-ZT13.NA.not_applicable-pval_05.csv')

gse_nash_plot(deq_nash_FANA_ZT13NAnomash)+ggtitle("ZT13 NA vs. FA non-MASH")+
  scale_y_discrete(position = "right") 
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.NA.kegg.pdf",height=6, width=12)

gse_nash_plot(deq_nash_FTNA_ZT13NAnomash)+ggtitle("ZT13 NA vs. FT non-MASH")+
  scale_y_discrete(position = "right") 
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT13.FT.Non_NASH-ZT13.NA.kegg.pdf",height=6, width=12)

gse_nash_plot(deq_nash_ZT113_NAnonmash)+ggtitle("ZT1 vs. ZT13 NA non-MASH")+
  scale_y_discrete(position = "right") 
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1.NA.not_applicable-ZT13.NA.not_applicable.kegg.pdf",height=6, width=12)

#looking at pathways enriched among all 3 pairwise condition comparisons
deq_nash_FANA_ZT13NAnomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.NA.not_applicable-pval_05.csv')
deq_nash_FTNA_ZT13NAnomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.Non_NASH-ZT13.NA.not_applicable-pval_05.csv')
deq_nash_FAFT_ZT13NAnomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.FT.Non_NASH-pval_05.csv')

allcond_nash<-gse_compareclust_allcond(deq_nash_FANA_ZT13NAnomash,deq_nash_FTNA_ZT13NAnomash,deq_nash_FAFT_ZT13NAnomash,
                                    "ZT13 non-MASH","ZT13 non-MASH","ZT13 non-MASH")
allcond_nash_df<-allcond_nash@compareClusterResult%>%
filter(!(category=="Human Diseases" | is.na(category)))
allcond_nash@compareClusterResult = allcond_nash@compareClusterResult[allcond_nash@compareClusterResult$ID %in% allcond_nash_df$ID, ]
dotplot(allcond_nash,x="direction", showCategory=200) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT13nonmash_allcond.compclust.nohumandisease.kegg.pdf",height=20, width=9)

#looking at pathways enriched among all 3 pairwise condition comparisons that are not shared by all
#also just NAFA and FAFT
deq_nash_FANA_ZT13NAnomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.NA.not_applicable-pval_05.csv')%>%
  filter(!(V1 %in% allcondpairwise))
  #filter(V1 %in% justNAFAFAFT)
deq_nash_FTNA_ZT13NAnomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.Non_NASH-ZT13.NA.not_applicable-pval_05.csv')%>%
  filter(!(V1 %in% allcondpairwise))
  #filter(V1 %in% justNAFAFAFT)
deq_nash_FAFT_ZT13NAnomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.FT.Non_NASH-pval_05.csv')%>%
  filter(!(V1 %in% allcondpairwise))
  #filter(V1 %in% justNAFAFAFT)

allcond_nash<-gse_compareclust_allcond(deq_nash_FANA_ZT13NAnomash,deq_nash_FTNA_ZT13NAnomash,deq_nash_FAFT_ZT13NAnomash,
                                       "ZT13 non-MASH","ZT13 non-MASH","ZT13 non-MASH")
allcond_nash<-gse_compareclust_allcond(deq_nash_FANA_ZT13NAnomash,deq_nash_FAFT_ZT13NAnomash,
                                       "ZT13 non-MASH","ZT13 non-MASH") #nothing enriched
allcond_nash_df<-allcond_nash@compareClusterResult%>%
  filter(!(category=="Human Diseases" | is.na(category)))
allcond_nash@compareClusterResult = allcond_nash@compareClusterResult[allcond_nash@compareClusterResult$ID %in% allcond_nash_df$ID, ]
dotplot(allcond_nash,x="direction", showCategory=30) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT13nonmash_allcond.compclust.nohumandisease.rmshared.kegg.pdf",height=15, width=9)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT13nonmash_allcond.compclust.nohumandisease.rmshared.top30.kegg.pdf",height=10, width=8)

deq_nash_mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT1.FA.NASH-ZT1.FT.NASH-pval_05.csv')
x<-gse_nash_plot(deq_nash_mash)
FAFT_nash_df<-x@result
# +ggtitle("ZT1 FA vs. FT MASH")+
#   scale_y_discrete(position = "right") 
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1.FA.NASH-ZT1.FT.NASH.kegg.pdf",height=6, width=12)

FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,"ZT1 MASH")
dotplot(FAFT_nash,x="direction", showCategory=100) +
  facet_grid(~condition) + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
FAFT_nash_df<-FAFT_nash@compareClusterResult


deq_nash_mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FT.NASH-pval_05.csv')
gse_nash_plot(deq_nash_mash)+ggtitle("ZT13 FA vs. FT MASH")+
  scale_y_discrete(position = "right")  #no pathways
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FT.NASH.kegg.pdf",height=6, width=12)

deq_nash_mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.FT.Non_NASH-pval_05.csv')
gse_nash_plot(deq_nash_mash)+ggtitle("ZT13 FA vs. FT non-MASH")+
  scale_y_discrete(position = "right")
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT13.FA.nonNASH-ZT13.FT.nonNASH.kegg.pdf",height=6, width=12)

deq_nash_ZT1mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT1.FA.NASH-ZT1.FT.NASH-pval_05.csv')
deq_nash_ZT13mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FT.NASH-pval_05.csv')
deq_nash_ZT13nonmash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.FT.Non_NASH-pval_05.csv')

FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,deq_nash_ZT13mash,
                            "ZT1 MASH","ZT13 non-MASH","ZT13 MASH")
# FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,
#                             "ZT1 MASH","ZT13 non-MASH")
dotplot(FAFT_nash,x="direction", showCategory=15) +
  facet_grid(~condition) + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1mashZT13nonmashFAFT.compclust_v2.pdf",height=10, width=9)

#plotting kegg
FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,deq_nash_ZT13mash,
                            "ZT1 MASH","ZT13 non-MASH","ZT13 MASH")
# FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,
#                             "ZT1 MASH","ZT13 non-MASH")
FAFT_nash_df<-FAFT_nash@compareClusterResult%>%
  filter(!(category=="Human Diseases" | is.na(category)))
FAFT_nash@compareClusterResult = FAFT_nash@compareClusterResult[FAFT_nash@compareClusterResult$ID %in% FAFT_nash_df$ID, ]
dotplot(FAFT_nash,x="direction", showCategory=15) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(angle = 45, vjust = 0.5))
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1mashZT13nonmashFAFT.compclust_v2.kegg.pdf",height=10, width=9)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1mashZT13nonmashFAFT.compclust_v2.nohumandisease.kegg.pdf",height=10, width=9)

#just overlap
deq_nash_ZT1mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT1.FA.NASH-ZT1.FT.NASH-pval_05.csv')%>%
  filter(V1 %in% just_shared)
deq_nash_ZT13mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FT.NASH-pval_05.csv')%>%
  filter(V1 %in% just_shared)
deq_nash_ZT13nonmash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.FT.Non_NASH-pval_05.csv')%>%
  filter(V1 %in% just_shared)

# FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,deq_nash_ZT13mash,
#                             "ZT1 MASH","ZT13 non-MASH","ZT13 MASH")
FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,
                            "ZT1 MASH","ZT13 non-MASH")
dotplot(FAFT_nash,x="direction", showCategory=50) +
  facet_grid(~condition) + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1mashZT13nonmashFAFT.compclust_v2.justshared.pdf",height=10, width=9)

#plotting kegg
# FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,deq_nash_ZT13mash,
#                             "ZT1 MASH","ZT13 non-MASH","ZT13 MASH")
FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,
                            "ZT1 MASH","ZT13 non-MASH")
FAFT_nash_df<-FAFT_nash@compareClusterResult%>%
  filter(!(category=="Human Diseases" | is.na(category)))
shared_pathway<-FAFT_nash_df$Description
FAFT_nash@compareClusterResult = FAFT_nash@compareClusterResult[FAFT_nash@compareClusterResult$ID %in% FAFT_nash_df$ID, ]
dotplot(FAFT_nash,x="direction", showCategory=20) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1mashZT13nonmashFAFT.compclust_v2.justshared.kegg.pdf",height=10, width=9)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1mashZT13nonmashFAFT.compclust_v2.justsharednohumandisease.kegg.pdf",height=6, width=9)


#just non-overlap
deq_nash_ZT1mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT1.FA.NASH-ZT1.FT.NASH-pval_05.csv')%>%
  filter(!(V1 %in% just_shared))
deq_nash_ZT13mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FT.NASH-pval_05.csv')%>%
  filter(!(V1 %in% just_shared))
deq_nash_ZT13nonmash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.FT.Non_NASH-pval_05.csv')%>%
  filter(!(V1 %in% just_shared))

# FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,deq_nash_ZT13mash,
#                             "ZT1 MASH","ZT13 non-MASH","ZT13 MASH")
FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,
                            "ZT1 MASH","ZT13 non-MASH")
dotplot(FAFT_nash,x="direction", showCategory=15) +
  facet_grid(~condition) + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1mashZT13nonmashFAFT.compclust_v2.nonoverlap.pdf",height=10, width=9)

#plotting kegg
# FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,deq_nash_ZT13mash,
#                             "ZT1 MASH","ZT13 non-MASH","ZT13 MASH")
# FAFT_nash<-gse_compareclust(deq_nash_ZT1mash,deq_nash_ZT13nonmash,
#                             "ZT1 MASH","ZT13 non-MASH")
FAFT_nash_df<-FAFT_nash@compareClusterResult%>%
  #filter(!(Description %in% shared_pathway))%>%
  filter(!(category=="Human Diseases" | is.na(category)))
FAFT_nash@compareClusterResult = FAFT_nash@compareClusterResult[FAFT_nash@compareClusterResult$ID %in% FAFT_nash_df$ID, ]
dotplot(FAFT_nash,x="direction", showCategory=15) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1mashZT13nonmashFAFT.compclust_v2.nonoverlap.kegg.pdf",height=10, width=9)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1mashZT13nonmashFAFT.compclust_v2.nonoverlap.rmhumandisease.kegg.pdf",height=8, width=9)
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1mashZT13nonmashFAFT.compclust_v2.nonoverlappathway.kegg.pdf",height=8, width=9)
#ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT1mashZT13nonmashFAFT.compclust_v2.nonoverlappathway.rmhumandisease.kegg.pdf",height=8, width=9)

###############################################################

#FAFT_nash_df<-FAFT_nash%>%as.data.frame()

nonNASH_bile<-allcond_nash_df%>%filter(condition=="ZT13 non-MASH" & grepl("Bile secretion",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="Bile secretion")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

nonNASH_insulin<-allcond_nash_df%>%filter(condition=="ZT13 non-MASH" & grepl("Insulin signaling pathway",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="Insulin signaling pathway")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

nonNASH_all<-FAFT_nash_df%>%filter(condition=="ZT13 non-MASH")%>% #& grepl("Insulin signaling pathway",Description)
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="Insulin signaling pathway")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

allcond_nash_df<-allcond_nash%>%as.data.frame()

nonNASHFT_insulin<-allcond_nash_df%>%filter(direction=="c-FT" & grepl("Insulin signaling pathway",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="Insulin signaling pathway")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

nonNASHFT_ampk<-allcond_nash_df%>%filter(direction=="c-FT" & grepl("AMPK",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="AMPK signaling pathway")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

nonNASHFT_mtor<-allcond_nash_df%>%filter(direction=="c-FT" & grepl("mTOR",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="mTOR signaling pathway")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

nonNASHFT_foxo<-allcond_nash_df%>%filter(direction=="c-FT" & grepl("FoxO",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="FoxO signaling pathway")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

md<-fread("transcriptomics/liver_metadata_wnashscore_fibste.tsv") %>%
  mutate(condition=ifelse(is.na(condition),"NA", condition),
         NASH_category= ifelse(NASH_category=="not_applicable", "Non_NASH",NASH_category))%>%
  mutate(mashZT=paste(timepoint,NASH_category,sep="_"))

#plot heatmap of these specific hits (Z-score)
genecounts<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.vsd.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  #filter(gene_id %in% nonNASHFT_ampk$gene_id)%>%
  #filter(gene_id %in% nonNASHFT_mtor$gene_id)%>%
  #filter(gene_id %in% nonNASHFT_foxo$gene_id)%>%
  filter(gene_id %in% nonNASH_insulin$gene_id)%>%
  gather(sampleid,TPM_counts,-gene_id) %>%
  left_join(.,md,by="sampleid") %>%
  #filter(mashZT=="ZT1_NASH"|mashZT=="ZT13_Non_NASH")%>%
  #filter(condition!="NA")%>%
  #left_join(.,sel_genes,by="gene_id")%>%
  group_by(gene_id,condition,mashZT)%>%dplyr::summarise(mn_TPM=mean(TPM_counts))%>%
  group_by(gene_id)%>%mutate(Zscore=(mn_TPM - mean(mn_TPM))/sd(mn_TPM))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         #mashZT=factor(mashZT,levels=c("ZT13_Non_NASH","ZT1_NASH")))%>%
         mashZT=factor(mashZT,levels=c("ZT1_Non_NASH","ZT13_Non_NASH","ZT1_NASH","ZT13_NASH")))%>%
  #GO=factor(GO,levels=c("lyase","ligase","ubiquitin","GTP")))%>% #"ubiquitin",
  #GO=factor(GO,levels=c("lyase","binding","ubiquitin")))%>% #"ubiquitin",
  separate(gene_id, c("gene_number","gene_symbol"), remove=FALSE)

# Run clustering
#genecounts_df<-fread("transcriptomics/data_files/liver_gene_tpm_matrix.csv")%>%
genecounts_df<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.vsd.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  filter(gene_id %in% nonNASH_insulin$gene_id)
  #filter(gene_id %in% nonNASHFT_ampk$gene_id)
  #filter(gene_id %in% nonNASHFT_mtor$gene_id)
  #filter(gene_id %in% nonNASHFT_foxo$gene_id)
genecounts_matrix<-genecounts_df%>%column_to_rownames("gene_id") %>%as.matrix()
genecounts_dendro <- as.dendrogram(hclust(d = dist(x = genecounts_matrix)))
dendro_plot <- ggdendrogram(data = genecounts_dendro, rotate = TRUE)

genecounts_order <- order.dendrogram(genecounts_dendro)
genecounts$gene_id <- factor(x = genecounts$gene_id,
                             levels = genecounts_df$gene_id[genecounts_order], 
                             ordered = TRUE)
#specify order
order_insulin<-c("Raf1","Prkar2a","Tsc1", "Rptor","Socs4",
                 "Shc1","Gsk3b","Phkb","Sos1","Mapk9","Akt2",
                 "Nras","Irs2","Ppp1cb","Phka2","Pygl","Pik3ca",
                 "Prkaa1","Mapk8","Braf","Sos2","Acaca","Insr",
                 "Cblb","Mknk1","Hk3","Mtor","Ppp1r3e","Rps6kb1","Prkacb",
                 "Pik3cb","Pde3b","Acacb","Rhoq","Prkcz","Gm5601","Pik3r1",
                 "Foxo1","Socs2","Pklr","Calml4","Gys2")
genecounts$gene_symbol<- factor(x = genecounts$gene_symbol,
                             levels = rev(order_insulin), 
                             ordered = TRUE)

plt<-ggplot(genecounts,aes(x=condition, y=gene_symbol)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+scale_y_discrete(position = "right") + 
  facet_grid(~mashZT,scales="free",space="free",switch = "y")+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),axis.text.y = element_text(size = 6),
        panel.spacing.y=unit(0.01, "lines"),strip.text.y = element_text(angle = 0))+
  scale_fill_viridis_c(option = "plasma")

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.heatmap.insulin_allcondZTmash.pdf", plot=plt,height=4, width=5)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.heatmap.insulin_ZT13nonmashZT1mash.pdf", plot=plt,height=4, width=3.5)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.heatmap.ampk_allcondZTmash.pdf", plot=plt,height=4, width=5)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.heatmap.mtor_allcondZTmash.pdf", plot=plt,height=4, width=5)
ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.heatmap.foxo_allcondZTmash.pdf", plot=plt,height=4, width=5)

###############################################################
deq_nash_ZT1mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT1.FA.NASH-ZT1.FT.NASH-pval_05.csv')
deq_nash_ZT13mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FT.NASH-pval_05.csv')
deq_nash_ZT13nonmash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.FT.Non_NASH-pval_05.csv')

list_venn <- list(ZT1_MASH = deq_nash_ZT1mash$V1,
                  ZT13_nonMASH = deq_nash_ZT13nonmash$V1,
                  ZT13_MASH = deq_nash_ZT13mash$V1)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

p<-venn.diagram(list_venn, height = 10,lty = 1,
                main="FA vs. FT",
                width = 10,lwd =1,filename = NULL)

grid.draw(p)
pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_ZT1ZT13nonmash.pdf")
grid.draw(p)
dev.off()

list_venn <- list(ZT1_MASH = deq_nash_ZT1mash$V1,
                  ZT13_nonMASH = deq_nash_ZT13nonmash$V1)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

just_shared<-all$`ZT1_MASH:ZT13_nonMASH`

p<-venn.diagram(list_venn, height = 10,lty = 1,
                main="FA vs. FT",
                width = 10,lwd =1,filename = NULL)

grid.draw(p)
pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_ZT1mashZT13nonmash.pdf")
grid.draw(p)
dev.off()


list_venn <- list(ZT13_MASH = deq_nash_ZT13mash$V1,
                  ZT13_nonMASH = deq_nash_ZT13nonmash$V1)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

p<-venn.diagram(list_venn, height = 10,lty = 1,
                main="FA vs. FT",
                width = 10,lwd =1,filename = NULL)

grid.draw(p)
pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_ZT13nonmashmash.pdf")
grid.draw(p)
dev.off()


deq_nash_FANA_ZT13NAnomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.NA.not_applicable-pval_05.csv')
deq_nash_FTNA_ZT13NAnomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.Non_NASH-ZT13.NA.not_applicable-pval_05.csv')
deq_nash_ZT13nonmash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.FT.Non_NASH-pval_05.csv')

list_venn <- list(NAFA = deq_nash_FANA_ZT13NAnomash$V1,
                  NAFT = deq_nash_FTNA_ZT13NAnomash$V1,
                  FAFT = deq_nash_ZT13nonmash$V1)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

p<-venn.diagram(list_venn, height = 10,lty = 1,
                main="ZT13 non-MASH",
                width = 10,lwd =1,filename = NULL)

justFAFT<-all$FAFT
justNAFAFAFT<-all$`NAFA:FAFT`
allcondpairwise<-all$`NAFA:NAFT:FAFT`

grid.draw(p)
pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_allconditions_ZT13nonmash.pdf")
grid.draw(p)
dev.off()

#upset plot
m = make_comb_mat(list_venn)
pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_allconditions_ZT13nonmash.upset.pdf",width=6,height=2.5)
UpSet(m, set_order=c("NAFA","NAFT","FAFT"))
dev.off()

list_venn <- list(NAFA = deq_nash_FANA_ZT13NAnomash$V1,
                  NAFT = deq_nash_FTNA_ZT13NAnomash$V1)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

p<-venn.diagram(list_venn, height = 10,lty = 1,
                main="ZT13 non-MASH",
                width = 10,lwd =1,filename = NULL)

grid.draw(p)
pdf(file="transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_ZT13nonmashnoFT.pdf")
grid.draw(p)
dev.off()

###############################################################

#volcano plot

#make volcano plots
annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")%>%
  filter(molecule_type=="protein_coding")

deq_nash_ZT1mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT1.FA.NASH-ZT1.FT.NASH.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%
  mutate(log2FoldChange=-1*log2FoldChange)%>%
  left_join(.,annot,by="stringtie_id")

p<-EnhancedVolcano(
  deq_nash_ZT1mash,
  lab = deq_nash_ZT1mash$gene_symbol,
  xlim=c(-10,10),
  ylim=c(0,40),
  title = 'FA vs. FT (ZT1 MASH)',
  subtitle=NA,
  caption=NA,
  x = "log2FoldChange",
  y = "padj",
  labSize = 0,
  pointSize = 1,
  colAlpha=1,
  col=c('gray37', 'gray37','blue','red'),
  pCutoff = 0.05,
  legendPosition = "none",
  border="full",
  gridlines.major=FALSE,
  gridlines.minor=FALSE
) +theme(text=element_text(size=7))

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2timepoint.condition.NASH_category-ZT1.FA.NASH-ZT1.FT.NASH.volplot_v3.pdf",height=5, width=5)

deq_nash_ZT13nonmash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.FT.Non_NASH.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%
  mutate(log2FoldChange=-1*log2FoldChange)%>%
  left_join(.,annot,by="stringtie_id")

p<-EnhancedVolcano(
  deq_nash_ZT13nonmash,
  lab = deq_nash_ZT13nonmash$gene_symbol,
   xlim=c(-10,10),
  ylim=c(0,40),
  title = 'FA vs. FT (ZT13 non-MASH)',
  subtitle=NA,
  caption=NA,
  x = "log2FoldChange",
  y = "padj",
  labSize = 0,
  colAlpha=1,
  pointSize = 1,
  col=c('gray37', 'gray37',"blue","red"),
  pCutoff = 0.05,
  legendPosition = "none",
  border="full",
  gridlines.major=FALSE,
  gridlines.minor=FALSE
) +theme(text=element_text(size=7))

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.FT.Non_NASH.volplot_v3.pdf",height=5, width=5)


deq_nash_ZT13mash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FT.NASH.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%
  mutate(log2FoldChange=-1*log2FoldChange)%>%
  left_join(.,annot,by="stringtie_id")

p<-EnhancedVolcano(
  deq_nash_ZT13mash,
  lab = deq_nash_ZT13mash$gene_symbol,
  xlim=c(-5,5),
  ylim=c(0,30),
  title = 'FA vs. FT (ZT13 MASH)',
  subtitle=NA,
  caption=NA,
  x = "log2FoldChange",
  y = "padj",
  labSize = 3,
  colAlpha=1,
  pointSize = 1,
  col=c('gray37', 'gray37',"blue","red"),
  pCutoff = 0.05,
  legendPosition = "none",
  border="full",
  gridlines.major=FALSE,
  gridlines.minor=FALSE
) +theme(text=element_text(size=7))

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT13.FA.NASH-ZT13.FT.NASH.volplot_v2.pdf",height=5, width=5)


deq_nash_NAFAZT13nomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.NA.not_applicable.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%
  mutate(log2FoldChange=log2FoldChange)%>%
  left_join(.,annot,by="stringtie_id")

p<-EnhancedVolcano(
  deq_nash_NAFAZT13nomash,
  lab = deq_nash_NAFAZT13nomash$gene_symbol,
  xlim=c(-10,10),
  ylim=c(0,40),
  title = 'NA vs. FA (ZT13 non-MASH)',
  subtitle=NA,
  caption=NA,
  x = "log2FoldChange",
  y = "padj",
  labSize = 0,
  colAlpha=1,
  pointSize = 1,
  col=c('gray37', 'gray37',"blue","red"),
  pCutoff = 0.05,
  legendPosition = "none",
  border="full",
  gridlines.major=FALSE,
  gridlines.minor=FALSE
) +theme(text=element_text(size=7))

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2.timepoint.condition.NASH_category-ZT13.FA.Non_NASH-ZT13.NA.not_applicable.volplot_v2.pdf",height=5, width=5)

deq_nash_NAFTZT13nomash<-fread('transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT13.FT.Non_NASH-ZT13.NA.not_applicable.csv')%>%
  dplyr::rename(stringtie_id=V1)%>%
  mutate(log2FoldChange=log2FoldChange)%>%
  left_join(.,annot,by="stringtie_id")

p<-EnhancedVolcano(
  deq_nash_NAFTZT13nomash,
  lab = deq_nash_NAFTZT13nomash$gene_symbol,
  xlim=c(-10,10),
  ylim=c(0,40),
  title = 'NA vs. FT (ZT13 non-MASH)',
  subtitle=NA,
  caption=NA,
  x = "log2FoldChange",
  y = "padj",
  labSize = 0,
  colAlpha=1,
  pointSize = 1,
  col=c('gray37', 'gray37',"blue","red"),
  pCutoff = 0.05,
  legendPosition = "none",
  border="full",
  gridlines.major=FALSE,
  gridlines.minor=FALSE
) +theme(text=element_text(size=7))

ggsave("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2..timepoint.condition.NASH_category-ZT13.FT.Non_NASH-ZT13.NA.not_applicable.volplot_v2.pdf",height=5, width=5)
