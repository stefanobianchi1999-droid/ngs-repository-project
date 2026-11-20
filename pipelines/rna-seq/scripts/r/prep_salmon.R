# Output: results/rna-seq/Exercise_GF/deseq2/three_groups_<tissue>/
#
# Run with the working directory set to the project root:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_exercise_gf.R")
# ============================================

# ============================================
# Exercise_GF: common setup (samples, tx2gene, Salmon import via tximport)
# Used by salmon_deseq2_exercise_gf.R (pairwise comparisons) and
# three_groups_exercise_gf.R (WT vs Bex3KO vs Bex3KO EE per tissue).
# Run with the working directory set to the project root.
# ============================================

library(tximport)
library(DESeq2)
library(ggplot2)

prep_salmon <- function(groups, tissues, reps, gtf_file, sample_dir) {
  # ================================================================
  # Step 1: sample table (from "Sample Information Collection Form")
  # ================================================================

  n_groups  <- length(groups)
  n_tissues <- length(tissues)
  combo     <- n_groups * n_tissues

  samples <- data.frame(
    np     = paste0("NP", 1:(combo * reps)),
    group  = rep(groups, each = n_tissues * reps),
    tissue = rep(tissues, each = reps, times = n_groups),
    rep    = rep(1:reps, times = combo)
  )

  samples$condition <- paste(samples$group, samples$tissue, sep = "_")
  samples$sample    <- paste(samples$condition, samples$rep, sep = "_")
  rownames(samples) <- samples$sample

  # ============================================
  # Step 2: transcript -> gene table (tx2gene) from GTF
  # ============================================
  # quote = "" keeps the attribute quotes (gene_id "ENSMUSG...";)

  gtf <- read.delim(
    gzfile(gtf_file),
    header = FALSE,
    comment.char = "#",
    quote = "",
    colClasses = c("character", "NULL", "character", rep("NULL", 5), "character")
  )
  colnames(gtf) <- c("chrom", "feature", "attributes")
  gtf <- gtf[gtf$feature == "transcript", ]

  # Extracts GTF attribute values (e.g. gene_id "ENSMUSG..."; -> ENSMUSG...)
  get_attr <- function(attrs, key) {
    sub(paste0(".*", key, ' "([^"]+)".*'), "\\1", attrs)
  }
  tx2gene <- data.frame(
    tx   = get_attr(gtf$attributes, "transcript_id"),
    gene = get_attr(gtf$attributes, "gene_id")
  )
  gene_names <- unique(data.frame(
    gene_id   = tx2gene$gene,
    chrom     = gtf$chrom,
    gene_name = ifelse(grepl("gene_name", gtf$attributes),
                        get_attr(gtf$attributes, "gene_name"), tx2gene$gene)
  ))
  symbol_of <- function(ids) gene_names$gene_name[match(ids, gene_names$gene_id)]
  chrom_of  <- function(ids) gene_names$chrom[match(ids, gene_names$gene_id)]

  # ============================================
  # Step 3: import all quant.sf via tximport (sums transcripts -> genes)
  # ============================================

  files <- file.path(sample_dir, samples$np, "quant.sf")
  names(files) <- samples$sample

  # quant.sf IDs are versioned (ENSMUST...1), the GTF's are not ->
  # ignoreTxVersion
  txi_all <- tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE
  )

  # Subset of txi keeping only some samples
  subset_txi <- function(txi, keep) {
    txi$abundance <- txi$abundance[, keep, drop = FALSE]
    txi$counts    <- txi$counts[, keep, drop = FALSE]
    txi$length    <- txi$length[, keep, drop = FALSE]
    txi
  }

  return(list(
    txi_all    = txi_all,
    samples    = samples,
    subset_txi = subset_txi,
    symbol_of  = symbol_of,
    chrom_of   = chrom_of,
    groups     = groups
  ))
}

# ============================================
# DESeq2 object for one tissue, built from the prep_salmon() output
# ============================================
dds_salmon_creation <- function(prep, tissue) {
  cd <- prep$samples[prep$samples$tissue == tissue, ]
  cd$group <- factor(cd$group, levels = prep$groups)
  txi <- prep$subset_txi(prep$txi_all, cd$sample)

  dds <- DESeqDataSetFromTximport(txi, colData = cd, design = ~ group)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  # I geni del cromosoma Y riflettono solo il sesso degli animali: li escludo
  dds <- dds[!prep$chrom_of(rownames(dds)) %in% "Y", ]

  # Test LRT: il gene cambia tra almeno due gruppi?
  DESeq(dds, test = "LRT", reduced = ~ 1)
}
