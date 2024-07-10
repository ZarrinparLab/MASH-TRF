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

##################################################################

md<-fread("transcriptomics/Panda_multiorgan_TRF/SraRunTable-liver-2.txt")%>%
  mutate(condition=ifelse(feeding_intervention=="ad-libitum feeding", "DIO","TRF"))
write.table(md,"transcriptomics/Panda_multiorgan_TRF/SraRunTable-liver-clean.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#convert data file
dat<-fread("transcriptomics/Panda_multiorgan_TRF/gene_count_matrix.csv")
write.table(dat,"transcriptomics/Panda_multiorgan_TRF/gene_count_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#filter for just protein coding sequence
annot<-fread("transcriptomics/data_files/stringtie_annotations.tsv")%>%
  filter(molecule_type=="protein_coding")
dat_sub<-dat%>%filter(gene_id %in% annot$stringtie_id)
write.table(dat_sub,"transcriptomics/Panda_multiorgan_TRF/gene_count_prot_matrix.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat_sub,"transcriptomics/Panda_multiorgan_TRF/gene_count_prot_matrix.csv",sep = ",",row.names = FALSE, quote=FALSE)

#Transcriptomics Reanalysis of STAM-HCC data (liver) --Deota 2023
##################################################################

#genelist<-c("Hoga1","Afmid","Shmt2", "Prodh", "Srm","Dhx58","Isg15","Ddx58","Ifih1","Rela","Rnf125")
genelist<-c("Acadl", "Acaa1a","Acox1","Ilk","Pparg","Slc27a1","Slc27a2","Slc27a5","Acox3")
results<-fread("transcriptomics/Deota2023_DEliver_results.csv")%>%
  filter(external_gene_name %in% genelist)%>%
  mutate(condition=ifelse(`logFC (TRF/ALF)`<0,"ALF","TRF"),
         external_gene_name = fct_reorder(external_gene_name, `logFC (TRF/ALF)`))
         # category=case_when(external_gene_name %in% c("Hoga1","Afmid","Shmt2") ~ "Glyoxylate and dicarboxylate metabolism",
         #                    external_gene_name %in% c("Srm", "Prodh") ~ "Arginine and proline metabolism",
         #                    external_gene_name %in% c("Dhx58","Isg15","Ddx58","Ifih1","Rela","Rnf125") ~ "RIG−I−like receptor signaling pathway"))%>%
  #mutate(category=factor(category,levels=c("Glyoxylate and dicarboxylate metabolism","Arginine and proline metabolism","RIG−I−like receptor signaling pathway")))

ggplot(results, aes(x =external_gene_name , y = `logFC (TRF/ALF)`, fill=condition)) + 
  geom_bar(stat="identity", alpha=0.7)+ theme_bw()+
  scale_fill_manual(values=c("#D55E00","#009E73"))+
  labs(title="Deota et al. 2023\nliver RNA ")+
  #facet_grid(category~., scale="free_y", space="free")+
  #scale_y_continuous(limits=c(-0.55,0.55))+
  theme(legend.position = "top")+
  coord_flip()
#ggsave("transcriptomics/Deota2023_DEresult_selgene_v2.pdf",height=4, width=4)
ggsave("transcriptomics/Deota2023_DEresult_selgene_v3.pdf",height=3, width=4)

genelist<-c("Acaa1a","Apoa1","Fabp2","Pck1","Ppara","Slc27a4","Acsl5","Acsl3","Acsl3")

results<-fread("transcriptomics/Deota2023_DEileum_results.csv")%>%
  filter(external_gene_name %in% genelist)%>%
  mutate(condition=ifelse(`logFC (TRF/ALF)`<0,"ALF","TRF"),
         external_gene_name = fct_reorder(external_gene_name, `logFC (TRF/ALF)`))

ggplot(results, aes(x =external_gene_name , y = `logFC (TRF/ALF)`, fill=condition)) + 
  geom_bar(stat="identity", alpha=0.7)+ theme_bw()+
  scale_fill_manual(values=c("#D55E00","#009E73"))+
  labs(title="Deota et al. 2023\nileum RNA ")+
  #facet_grid(category~., scale="free_y", space="free")+
  #scale_y_continuous(limits=c(-0.55,0.55))+
  theme(legend.position = "top")+
  coord_flip()
ggsave("transcriptomics/Deota2023_DEresultileum_selgene_v2.pdf",height=4, width=4)

##################################################################
#just deota liver results 

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

gse_compareclust_ZT<-function(df1,cond1){
  gene_list1<-get_genelist(df1)
  
  mydf1 <- data.frame(Entrez=names(gene_list1), FC=gene_list1)%>%
    mutate(direction=ifelse(FC<0,"a-ZT13_non-MASH","b-ZT1-MASH"),
           condition=cond1)
  
  formula_res <- compareCluster(Entrez~direction+condition, organism = "mmu", data=mydf1,
                                fun="enrichKEGG",pvalueCutoff=0.05,pAdjustMethod="fdr")
  return(formula_res)
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
  
  p<-dotplot(kk2, showCategory=30, split=".sign") + facet_grid(.~.sign)
  return(p)
}


##################################################################

#MASH TRF 
MASHTRF_results<-fread("transcriptomics/NASH_liver_deseq2/other_comparisons/NASH_liver_deseq2_.lfchange.timepoint.condition.NASH_category-ZT1.FT.NASH-ZT13.FT.Non_NASH-pval_05.csv")%>%
  separate(V1,c("geneid","gene"),sep="\\|",remove=FALSE)

results<-fread("transcriptomics/Panda_multiorgan_TRF/Deota_liver_deseq2_.lfchange.timepoint.condition-ZT26.TRF-ZT38.TRF-pval_05.csv")%>%
  separate(V1,c("geneid","gene"),sep="\\|",remove=FALSE)

list_venn <- list(mystudy = MASHTRF_results$gene,
                  deota = results$gene)
ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

MASHTRF_overlap<-MASHTRF_results%>%filter(gene %in% all$mystudy)

ZT_nash<-gse_compareclust_ZT(MASHTRF_overlap,"ZT13 non-MASH vs. ZT1 MASH")
dotplot(ZT_nash,x="direction", showCategory=30) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
gse_nash_plot(MASHTRF_overlap)+labs(title="ZT1MASHvZT13NonMASH (subset ZT26vZT38 TRF)")

ggsave("transcriptomics/Panda_multiorgan_TRF/keggenrich_ZT1MASHZT13NonMASH_alsoinDeota_kegg.pdf",height=4, width=8)
ggsave("transcriptomics/Panda_multiorgan_TRF/keggenrich_ZT1MASHZT13NonMASH_justinourStudy_kegg.pdf",height=10, width=8)

#plot overlap

p<-venn.diagram(list_venn, fill = c("#D55E00","#009E73"),height = 10,lty = 0, 
                main="ZT1MASHvZT13NonMASH v. ZT26vZT38 TRF",
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="transcriptomics/Panda_multiorgan_TRF/venn_TRFs_ZT1MASHZT13NonMASHvZT26ZT38TRF.pdf")
grid.draw(p)
dev.off()


#looking at insulin genes
insulin_genes<-c("Raf1","Prkar2a","Tsc1", "Rptor","Socs4",
                 "Shc1","Gsk3b","Phkb","Sos1","Mapk9","Akt2",
                 "Nras","Irs2","Ppp1cb","Phka2","Pygl","Pik3ca",
                 "Prkaa1","Mapk8","Braf","Sos2","Acaca","Insr",
                 "Cblb","Mknk1","Hk3","Mtor","Ppp1r3e","Rps6kb1","Prkacb",
                 "Pik3cb","Pde3b","Acacb","Rhoq","Prkcz","Gm5601","Pik3r1",
                 "Foxo1","Socs2","Pklr","Calml4","Gys2")

results<-fread("transcriptomics/Panda_multiorgan_TRF/Deota_liver_deseq2_.lfchange.timepoint.condition-ZT26.TRF-ZT38.TRF.csv")%>%
  separate(V1,c("geneid","gene"),sep="\\|",remove=FALSE)%>%
  filter(gene %in% insulin_genes)

