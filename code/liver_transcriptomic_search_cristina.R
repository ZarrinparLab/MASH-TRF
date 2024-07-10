setwd("/mnt/zarrinpar/scratch/sfloresr/NASH_KF")

library(tidyverse)
library(data.table)

###############################################################
M4_pathway<-c("Camk2b","Chrm4","Itpr1","Cacna1a","Pik3cd","Cacna1c","Pik3r1","Adcy1","Adcy7","Pik3r6","Pik3cg","Pik3r5","Adcy5","Gnai2","Gngt2","Gng2","Gng4","Gng7","Akt3","Fyn","
              Kcnj3","Chrnb2","Kcnj6","Prkcb","Gnao1","Slc5a7","Camk4","Bcl2","Kcnq5","Gnb5","Plcb2")

Pi3k_Akt<-c("Csf3r","Flt3","Flt4","Itgb3","Tnc","Pik3cd","Pik3cg","Ghr","Gngt2","Fgf9","Tnn","Myc","Akt3","Myb","Tnr","Itgb7","Jak3","Itga4","Magi2","Rps6","Itga1","Ppp2r5c","Pgf","Col4a4","
Col6a2","Kit","Col6a5","Itga5","Itga9","Csf1r","Tnxb","Pik3r1","Pik3r6","Pik3r5","Vtn","Reln","Gng2","Gng4","Cd19","Pdgfd","Gng7","Ntrk2","Angpt1","Vegfc","Vegfd","Igf1","
Kitlg","Col9a1","Bcl2","Col9a2","Pkn1","Gnb5","Fgfr4","Il7r","Pik3ap1","Fgfr1")

Actin_Cytoskeleton<-c("Cyfip2","Chrm4","Itgb3","Itgb2","Cxcr4","Pik3cd","Itgal","Mylk","Slc9a1","Fgd3","Cfl2","Cfl1","Itgax","Pip4k2a","Myh14","Nckap1l","Myh11","Itgb7","Pip4k2c","Rac1","Itga1","
Cxcl12","Ppp1r12b","Arhgef6","Fgfr1","Itga9")

Dendritic_Cells<-c("Abca6","Abcc3","Adam8","Adgre1","Adgre2","Adgrg3","Adgrg5","Aldh2","Alox5","Anpep","Aoah","Apoc1","Arhgap30","Bach2","Bin2","Btk","C16orf54","C3","C3ar1","
Cadm1","Card11","Cass4","Ccdc141","Ccl18","Ccl21","Ccr5","Cd101","Cd19","Cd200r1","Cd22","Cd226","Cd300lf","Cd37","Cd38","Cd48","Cd52","Cd84","Cd96","Clec10a","
Clec12a","Clec2d","Clec4f","Clec7a","Clnk","Col19a1","Col6a5","Colec12","Cpa3","Cped1","Cpvl","Cr1","Cx3cr1","Cxcl12","Cyth4","Cytip","Dock2","Dock8","Dynlt1","
Entpd1","Evi2a","Evi2b","Fbxo6","Fermt3","Fgd2","Fgd3","Fgr","Fxyd6","Gab3","Gapt","Gfra1","Gng7","Gpr132","Gpr171","Gpr65","Gzmk","Hck","Hcls1","Hdc","Hhex","Hla-Doa","Hvcn1","Ido1","Igsf6","Ikzf1","Ikzf3","Il16","Il27ra","Ipcef1","Itgal","Itgax","Itgb2","Itgb7","Itk","Jak3","Jaml","Klrb1","Kxd1","Lair1","Lat2","Lcp2","Lgals2","Lilrb2","
Lrmp","Lrrk2","Lsp1","Ly9","Macrod2","Map4k1","Milr1","Mmp9","Ms4a2","Msr1","Myo1f","Nampt","Nckap1l","Nfam1","Nlrc3","P2rx1","P2rx7","P2ry13","P2ry14","
Parp15","Parvg","Pde1b","Pde3b","Pgm5","Pik3cg","Pik3r5","Pik3r6","Pip4k2a","Plcb2","Pld4","Ppp1r1a","Pram1","Prf1","Pstpip2","Ptafr","Ptgdr","Ptgds","Ptprc","
Rasgrp1","Rasgrp3","Rcsd1","Renbp","Rgs18","Rtn1","S1pr4","Samsn1","Scimp","Siglec10","Siglec11","Siglec6","Sla","Slc24a4","Slco5a1","Smco4","Snx20","Tbx21","
Tlr10","Tlr6","Tmc8","Tnfrsf13c","Tns1","Tpsab1","Traf3ip3","Uevld","Vipr2","Wdfy4","Xcr1","Znf366","Znf831,
Irf4"," Cd209"," Nrp1"," Batf3"," Zbtb46"," Irf8"," Btla"," Cadm1"," Cd8a"," Itgae"," Itgax"," Xcr1"," Irf4"," Cd163"," Clec10a"," Notch2"," Itgam"," Sirpa"," H2"," Cx3cr1"," 
Cd1c"," Il23"," Cxcl13"," Il12"," Fcer1a"," Mrc1"," Itggb7"," Ccr9"," Flt3l"," Tlr5"," Ptprc"," Ccr5"," Cxcr4"," Sdc3"," Siglec10"," Siglec11"," Cd209")

Amp<-c("Defa1"," Defa2"," Defa3"," Defa4"," Defa5"," Defa6"," Defb1-4"," Camp"," Cramp"," Lyz"," Pla2g2a"," Ang"," Bpi"," Ccl20"," Retnlb"," Tlr10"," Trl6"," Tlr"," Myd88"," Nod1"," Nod2"," Birc2"," Birc3"," Nlrc4"," Nlrc3")

###############################################################
#liver
md<-fread("transcriptomics/liver_metadata_wnashscore_fibste_wmashZT.tsv")%>%
  mutate(condition=ifelse(is.na(condition), "NA", condition))%>%
  filter(condition!="FT" & timepoint=="ZT1")
write.table(md,"transcriptomics/files_for_cristina/liver/metadata_liver_ZT1_FANA.txt",sep = "\t",row.names = FALSE,quote=FALSE)

dsq2<-fread("transcriptomics/NASH_liver_deseq2/NASH_liver_deseq2_.lfchange.condition.timepoint-FA.ZT1-NA.ZT1.csv") %>% # ZT1_FA_MASH/ZT1_NA_nonMASH
  separate(V1,c("gene_id","gene_name"))%>%
  filter(padj<0.05)
write.table(dsq2,"transcriptomics/files_for_cristina/liver/deseq_liver_results_ZT1_FA_MASHoverZT1_NA_NonMASH.txt",sep = "\t",row.names = FALSE,quote=FALSE)

M4_pathway_hits<-dsq2%>%filter(gene_name %in% M4_pathway)%>% #28/31 present
  filter(padj<0.05) #6 sig diff
write.table(M4_pathway_hits,"transcriptomics/files_for_cristina/liver/M4_pathway_liver_results.txt",sep = "\t",row.names = FALSE,quote=FALSE)

Pi3k_Akt_hits<-dsq2%>%filter(gene_name %in% Pi3k_Akt)%>% #55/56 present
  filter(padj<0.05) #9 sig diff
write.table(Pi3k_Akt_hits,"transcriptomics/files_for_cristina/liver/Pi3k_Akt_liver_results.txt",sep = "\t",row.names = FALSE,quote=FALSE)

Actin_Cytoskeleton_hits<-dsq2%>%filter(gene_name %in% Actin_Cytoskeleton)%>% #26/26 present
  filter(padj<0.05) #7 sig diff
write.table(Actin_Cytoskeleton_hits,"transcriptomics/files_for_cristina/liver/Actin_Cytoskeleton_liver_results.txt",sep = "\t",row.names = FALSE,quote=FALSE)

Dendritic_Cells_hits<-dsq2%>%filter(gene_name %in% Dendritic_Cells)%>% #142/207 present
  filter(padj<0.05) #23 sig diff
write.table(Dendritic_Cells_hits,"transcriptomics/files_for_cristina/liver/Dendritic_Cells_liver_results.txt",sep = "\t",row.names = FALSE,quote=FALSE)

Amp_hits<-dsq2%>%filter(gene_name %in% Amp)%>% #0/25 present
  filter(padj<0.05) #0 sig diff
write.table(Amp_hits,"transcriptomics/files_for_cristina/liver/Amp_liver_results.txt",sep = "\t",row.names = FALSE,quote=FALSE)

###############################################################
#ileum

md<-fread("transcriptomics/ileum_md_with_nash_score.tsv")%>%
  mutate(condition=ifelse(is.na(condition), "NA", condition),
         NASH_category=ifelse(NASH_category=="not_applicable","Non_NASH","NASH"))%>%
  filter(condition!="FT" & timepoint=="ZT1")
write.table(md,"transcriptomics/files_for_cristina/ileum/metadata_ileum_ZT1_FANA.txt",sep = "\t",row.names = FALSE,quote=FALSE)

dsq2<-fread("transcriptomics/NASH_ileum_deseq2/lowthreshold3_results/NASH_ileum_deseq2_.lfchange.condition.timepoint-FA.ZT1-NA.ZT1.csv") %>% # ZT1_FA_MASH/ZT1_NA_nonMASH
  separate(V1,c("gene_id","gene_name"))%>%
  filter(padj<0.05)
write.table(dsq2,"transcriptomics/files_for_cristina/ileum/deseq_ileum_results_ZT1_FA_MASHoverZT1_NA_NonMASH.txt",sep = "\t",row.names = FALSE,quote=FALSE)

M4_pathway_hits<-dsq2%>%filter(gene_name %in% M4_pathway)%>% #28/31 present
  filter(padj<0.05) #0 sig diff
write.table(M4_pathway_hits,"transcriptomics/files_for_cristina/ileum/M4_pathway_ileum_results.txt",sep = "\t",row.names = FALSE,quote=FALSE)

Pi3k_Akt_hits<-dsq2%>%filter(gene_name %in% Pi3k_Akt)%>% #51/56 present
  filter(padj<0.05) #4 sig diff
write.table(Pi3k_Akt_hits,"transcriptomics/files_for_cristina/ileum/Pi3k_Akt_ileum_results.txt",sep = "\t",row.names = FALSE,quote=FALSE)

Actin_Cytoskeleton_hits<-dsq2%>%filter(gene_name %in% Actin_Cytoskeleton)%>% #26/26 present
  filter(padj<0.05) #0 sig diff
write.table(Actin_Cytoskeleton_hits,"transcriptomics/files_for_cristina/ileum/Actin_Cytoskeleton_ileum_results.txt",sep = "\t",row.names = FALSE,quote=FALSE)

Dendritic_Cells_hits<-dsq2%>%filter(gene_name %in% Dendritic_Cells)%>% #134/207 present
  filter(padj<0.05) #19 sig diff
write.table(Dendritic_Cells_hits,"transcriptomics/files_for_cristina/ileum/Dendritic_Cells_ileum_results.txt",sep = "\t",row.names = FALSE,quote=FALSE)

Amp_hits<-dsq2%>%filter(gene_name %in% Amp)%>% #0/25 present
  filter(padj<0.05) #0 sig diff
write.table(Amp_hits,"transcriptomics/files_for_cristina/ileum/Amp_ileum_results.txt",sep = "\t",row.names = FALSE,quote=FALSE)


