#!/usr/bin/env Rscript
options(warn = -1)
library(argparser, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(gtools, quietly = TRUE)
suppressPackageStartupMessages(library(DESeq2, quietly = TRUE))
library(tximport, quietly = TRUE)
library(tximportData, quietly = TRUE)
library(BiocParallel)
library(EnhancedVolcano)
options(warn = 0)


main <- function(script) {
  opts <- getargs()

  # opts <- getargs(
  #  strsplit(
  #    '-m md_with_nash_score.tsv -c NASH_category -f hisat2/ileum.nash.pruned_min3_conditions_w_3_samples_gt_100.gene_count_matrix.csv -o hisat2/ileum_nash_pruned_min3_3_samples.deseq2/ileum.nash.pruned --mincount 100',
  #    ' ',
  #   )
  #  )
  register(MulticoreParam(opts$threads))

  # file names for output
  out_lf_change_tmpl <- sprintf("%s.lfchange.%%s.csv", opts$output)
  out_volcano_tmpl <- sprintf("%s.lfchange.%%s.svg", opts$output)
  out_counts <- sprintf("%s.counts.csv", opts$output)
  out_ntd <- sprintf("%s.ntd.csv", opts$output)
  out_vsd <- sprintf("%s.vsd.csv", opts$output)
  out_rld <- sprintf("%s.rld.csv", opts$output)
  out_dds <- sprintf("%s.dds.rds", opts$output)
  out_md <- sprintf("%s.metadata.csv", opts$output)
  out_pca <- sprintf("%s.pca.pdf", opts$output) 
  out_pca_tb <- sprintf("%s.pca.csv", opts$output)  

  # load count table
  tbl <- read.csv(opts$file, header = TRUE, row.names = 1, check.names = FALSE)
  # read metadata and sort columns by import order
  md_path <- opts$metadata
  md <- read.table(md_path, header = TRUE, sep = "\t", row.names = 1, na.strings = c(""), check.names = FALSE)

  # drop count table columns not in metadata
  tbl <- tbl[,intersect(rownames(md), colnames(tbl))]
  print(colnames(tbl))
  # reorder/drop metadata not in count tabl
  md <- md[colnames(tbl), ] # reorder mdat based on counts
  print(rownames(md))
    
   # change to vectors
  columns_compare <- as.vector(unlist(strsplit(opts$columns, ",")))
  md[columns_compare] <- lapply(md[columns_compare], factor)

  # create new column with all experimental variables
  if (length(columns_compare) > 1) {
    new_col_name <- paste(columns_compare, sep = ".", collapse = ".")
    md <- unite(md, {{ new_col_name }}, all_of(columns_compare), sep = ".", remove = FALSE)
    md[, new_col_name] <- as.factor(md[, new_col_name])
  } else {
    new_col_name <- columns_compare[1]
  }

  # set the condition formula
  condition <- as.formula(paste("~", new_col_name, sep = ""))

  # load the dds
  dds <- DESeqDataSetFromMatrix(countData = tbl, colData = md, design = condition)

  # filter out <5 counts or <5 samples
  dds <- estimateSizeFactors(dds)
  idx <- rowSums(counts(dds, normalized = TRUE) >= opts$mincount) >= 5
  dds <- dds[idx, ]

  # run Deseq2
  dds <- DESeq(dds, parallel = TRUE)

  # load and add annotations to dds
  annotations <- read.table(opts$annotations, sep = "\t", header = TRUE, row.names = 1, na.strings = c(""))

  # save accessory tables
  cnt <- counts(dds, normalized = TRUE)
  # cnt <- merge(
  #  as.data.frame(cnt), annotations, by=0, all.x=TRUE
  # ) %>% column_to_rownames('Row.names')
  write.csv(cnt, file = out_counts)

  vsd <- vst(dds, blind = FALSE)
  # vsd_a <- merge(
  #  as.data.frame(
  #    assay(vsd)), annotations, by=0, all.x=TRUE
  #  ) %>% column_to_rownames('Row.names')
  write.csv(assay(vsd), file = out_vsd)

  ntd <- normTransform(dds)
  # ntd_a <- merge(
  #  as.data.frame(assay(ntd)), annotations, by=0, all.x=TRUE
  #  ) %>% column_to_rownames('Row.names')
  write.csv(assay(ntd), file = out_ntd)

  if (opts$rld) {
    rld <- rlog(dds, blind = FALSE)
    # rld_a = merge(as.data.frame(assay(rld)), annotations, by=0, all.x=TRUE
    # ) %>% column_to_rownames('Row.names')
    write.csv(assay(rld), file = out_rld)
  }

  # write out pairwise comparisons
  sample_types <- as.character(unique(md[, new_col_name]))
  pairwise_samples <- combinations(length(sample_types), 2, sample_types)
  for (n in 1:nrow(pairwise_samples)) {
    v <- prepend(pairwise_samples[n, ], {{ new_col_name }})
    res <- results(dds, contrast = v, parallel = TRUE)
    res_df <- as.data.frame(res)
    # res_df <- res_df[!is.na(res_df$padj), ]
    # res_df <- merge(
    #  res_df, annotations, by=0, all.x=TRUE
    #  ) %>% column_to_rownames('Row.names')
    res_df_sig <- subset(res_df, padj <= 0.05)
    print(c(v, nrow(res_df_sig)))
    write.csv(res_df, sprintf(out_lf_change_tmpl, paste(v, collapse = "-")))
    write.csv(res_df_sig, sprintf(out_lf_change_tmpl, paste(c(v, "pval_05"), collapse = "-")))
    svg(sprintf(out_volcano_tmpl, paste(v, collapse = "-")))
    print(EnhancedVolcano(
      res_df,
      lab = rownames(res_df),
      x = "log2FoldChange",
      y = "padj",
      pCutoff = 0.05
    ))
    dev.off()
  }
}


getargs <- function(argv = commandArgs(trailingOnly = TRUE)) {
  p <- arg_parser("Run DESeq Basic") ## TODO program description
  p <- add_argument(p, "--metadata", "-m", help = "metadata file")
  p <- add_argument(p, "--columns", "-c", help = "comma(,) delimited metadata columns for comparison")
  p <- add_argument(p, "--mincount", default = 5, help = "minimum count per sample to include")
  p <- add_argument(p, "--file", "-f", help = "counts table (csv)")
  p <- add_argument(p, "--output", "-o", help = "output prefix")
  p <- add_argument(p, "--threads", "-j", default = 24, help = "threads for DESeq")
  p <- add_argument(p, "--annotations", "-a",
    default = "/mnt/zarrinpar/Pynchon/Databases/mouse_genome/gene_annotations.tsv",
    help = "gene annotation table (gene_id, gene_name, description)"
  )
  p <- add_argument(p, "--rld", flag = TRUE, help = "Create RLD file (slow)")
  opts <- parse_args(p, argv)
  return(opts)
}

# ----------- Helper methods -------------

thisfile <- function() { # find the path to the source file
  # returns path to the script
  if (!is.null(res <- thisfile_source())) {
    res
  } else if (!is.null(res <- thisfile_rscript())) {
    res
  } else if (!is.null(res <- thisfile_knit())) {
    res
  } else if (!is.null(res <- thisfile_rstudio_source())) {
    res
  } else if (!is.null(res <- thisfile_rstudio_run())) {
    res
  } else {
    NULL
  }
}

thisfile_source <- function() { # path if program is loaded by source() in another R script
  for (i in -(1:sys.nframe())) {
    if (identical(args(sys.function(i)), args(base::source))) {
      return(normalizePath(sys.frame(i)$ofile))
    }
  }
  NULL
}

thisfile_rstudio_source <- function() { # path if program is loaded by source() in RStudio
  for (i in -(1:sys.nframe())) {
    if (identical(args(sys.function(i)), args(base::source))) {
      return(normalizePath(sys.frame(i)$fileName))
    }
  }
  NULL
}

thisfile_rscript <- function() { # path if run as script from command line
  cmd_args <- commandArgs(trailingOnly = FALSE)
  cmd_args_trailing <- commandArgs(trailingOnly = TRUE)
  leading_idx <-
    seq.int(from = 1, length.out = length(cmd_args) - length(cmd_args_trailing))
  cmd_args <- cmd_args[leading_idx]
  res <- gsub("^(?:--file=(.*)|.*)$", "\\1", cmd_args)
  # If multiple --file arguments are given, R uses the last one
  res <- tail(res[res != ""], 1)
  if (length(res) > 0) {
    return(res)
  }
  NULL
}

thisfile_knit <- function() { # path is loaded via knitr
  if (requireNamespace("knitr")) {
    return(knitr::current_input())
  }
  NULL
}

thisfile_rstudio_run <- function() { # path if run via RStudio run()
  if (requireNamespace("rstudioapi")) {
    return(normalizePath(rstudioapi::getActiveDocumentContext()$path))
  }
  NULL
}
script.path <- thisfile()
script.dir <- dirname(script.path)

main(script)
