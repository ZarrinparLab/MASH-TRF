library(DESeq2)

coldata<-fread("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/transcriptomics/liver_metadata_wnashscore_fibste.tsv")%>%
  mutate(NASH_category=factor(NASH_category,levels=c("Non_NASH","NASH","Not_applicable")))
coldata_FA<-coldata%>%filter(condition=="FA")
coldata_FT<-coldata%>%filter(condition=="FT")
cts<-fread("/mnt/zarrinpar/scratch/sfloresr/NASH_KF/transcriptomics/data_files/liver_gene_count_prot_matrix.csv")
cts_FA<-cts%>%dplyr::select(gene_id,coldata_FA$sampleid)%>%column_to_rownames("gene_id")
cts_FT<-cts%>%dplyr::select(gene_id,coldata_FT$sampleid)%>%column_to_rownames("gene_id")
  


ddsFA <- DESeqDataSetFromMatrix(countData = cts_FA,
                              colData = coldata_FA,
                              design= ~ NASH_category)
ddsFA <- DESeq(ddsFA)
resultsNames(ddsFA) # lists the coefficients
res <- results(ddsFA)

df<-as.data.frame(res)%>%rownames_to_column("gene_id")
