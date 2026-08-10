#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: run_gcm1_deseq2.R INPUT_DIRECTORY OUTPUT_DIRECTORY")
}

if (getRversion() != package_version("4.2.2")) {
  stop(sprintf("exact R 4.2.2 required; found %s", getRversion()))
}

suppressPackageStartupMessages(library(DESeq2))
if (packageVersion("DESeq2") != package_version("1.36.0")) {
  stop(sprintf("exact DESeq2 1.36.0 required; found %s", packageVersion("DESeq2")))
}

input_directory <- args[[1]]
output_directory <- args[[2]]
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

for (state in c("EVT", "ST")) {
  counts <- read.delim(
    file.path(input_directory, paste0(state, ".counts.tsv")),
    row.names = 1,
    check.names = FALSE
  )
  metadata <- read.delim(
    file.path(input_directory, paste0(state, ".metadata.tsv")),
    row.names = 1,
    check.names = FALSE
  )
  if (!identical(colnames(counts), rownames(metadata))) {
    stop(paste(state, "count columns and metadata rows differ"))
  }
  metadata$condition <- relevel(factor(metadata$condition), ref = "WT")
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = metadata,
    design = ~condition
  )
  keep <- rowSums(counts(dds)) > 0
  dds <- DESeq(dds[keep, ], quiet = FALSE)
  result <- results(dds, contrast = c("condition", "GCM1_KO", "WT"))
  table <- data.frame(gene = rownames(result), as.data.frame(result), check.names = FALSE)
  write.table(
    table,
    file.path(output_directory, paste0(state, ".deseq2.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

session_file <- file.path(output_directory, "sessionInfo.txt")
sink(session_file)
sessionInfo()
sink()
