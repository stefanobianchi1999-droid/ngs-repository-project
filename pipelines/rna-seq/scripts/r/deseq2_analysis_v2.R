# Output: results/rna-seq/Exercise_GF/deseq2/three_groups_<tessuto>/
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_exercise_gf.R")
# ============================================

# ============================================
# Exercise_GF: setup comune (campioni, tx2gene, import Salmon con tximport)
# Usato da salmon_deseq2_exercise_gf.R (confronti a coppie) e
# three_groups_exercise_gf.R (WT vs Bex3KO vs Bex3KO EE per tessuto).
# Da lanciare con la working directory nella root del progetto.
# ============================================
library(tximport)
library(DESeq2)
library(ggplot2)

salmon_dir <- "results/rna-seq/Exercise_GF/salmon"

deseq_dir  <- "results/rna-seq/Exercise_GF/deseq2"


group_levels <- snakemake@params[["groups"]]
padj_cut     <- snakemake@params[["padj_cut"]]
tissue       <- snakemake@wildcards[["tissue"]]

# =============================================================================
# Step 1: sample table (from "Sample Information Collection Form")
# =============================================================================

groups <- snakemake@params[["groups"]]
tissues <- snakemake@params[["tissues"]]
reps <- snakemake@params[["reps"]]

n_groups <- length(groups)
n_tissues <- length(tissues) 
combo <- n_groups * n_tissues

samples <- data.frame(
  np     = paste0("NP", 1:(combo * reps)),
  group  = rep(groups, each = n_tissues * reps),
  tissue = rep(tissues,  each = reps, times = n_groups),
  rep    = rep(1:reps, times = combo)
)

samples$condition <- paste(samples$group, samples$tissue, sep = "_")
samples$sample <- paste(samples$condition, samples$rep, sep = "_")
rownames(samples) <- samples$sample

# ============================================
# Step 2: Trascriptome table -> gene (tx2gene) from GTF
# ============================================
# quote = "" tiene le virgolette degli attributi (gene_id "ENSMUSG...";)

gtf_file   <- snakemake@params[["gtf_file"]]

gtf <- read.delim(
  gzfile(gtf_file),
  header = FALSE,
  comment.char = "#",
  quote = "",
  colClasses = c("character", "NULL", "character", rep("NULL", 5), "character")
)
colnames(gtf) <- c("chrom", "feature", "attributes")
gtf <- gtf[gtf$feature == "transcript", ]

# This extracts GTF attribute values (es. gene_id "ENSMUSG..."; -> ENSMUSG...)
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
# Step 3: import of all quant.sf using tximport (it sums transcripts -> genes)
# ============================================
files <- file.path(salmon_dir, samples$np, "quant.sf")
names(files) <- samples$sample
# quant.sf ha gli ID con versione (ENSMUST...1), il GTF senza -> ignoreTxVersion
txi_all <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)

# Sottoinsieme di txi con solo alcuni campioni
subset_txi <- function(txi, keep) {
  txi$abundance <- txi$abundance[, keep, drop = FALSE]
  txi$counts    <- txi$counts[, keep, drop = FALSE]
  txi$length    <- txi$length[, keep, drop = FALSE]
  txi
}

group_levels <- c("WT_Control", "Bex3KO", "Bex3KO_EE")
group_colors <- c(WT_Control = "grey40", Bex3KO = "firebrick", Bex3KO_EE = "steelblue")
padj_cut <- 0.05

# ============================================
# Classificazione dei pattern a partire dai tre contrasti a coppie
#   ko = Bex3KO vs WT, ee_wt = Bex3KO EE vs WT, ee_ko = Bex3KO EE vs Bex3KO
# ============================================
classify_pattern <- function(ko, ee_wt, ee_ko) {
  sig <- function(x) !is.na(x$padj) & x$padj < padj_cut
  s_ko <- sig(ko); s_eewt <- sig(ee_wt); s_eeko <- sig(ee_ko)
  d_ko <- sign(ko$log2FoldChange); d_eewt <- sign(ee_wt$log2FoldChange); d_eeko <- sign(ee_ko$log2FoldChange)

  p <- rep("Other", length(s_ko))
  # alterato nel KO, l'EE lo riporta verso il WT
  p[s_ko & s_eeko & d_eeko == -d_ko & !s_eewt]                    <- "KO altered, rescued by EE"
  p[s_ko & s_eeko & d_eeko == -d_ko & s_eewt & d_eewt == d_ko]     <- "KO altered, partially rescued"
  # alterato nel KO, l'EE non lo corregge
  p[s_ko & !s_eeko & s_eewt & d_eewt == d_ko]                     <- "KO altered, not rescued"
  # alterato nel KO, l'EE sta a meta' (nessun contrasto con l'EE significativo)
  p[s_ko & !s_eeko & !s_eewt]                                     <- "KO altered, EE intermediate"
  # non alterato nel KO, cambiato solo dall'EE
  p[!s_ko & s_eeko & s_eewt & d_eeko == d_eewt]                   <- "EE-specific"
  # KO ed EE si discostano dal WT in direzioni opposte, nessuno dei due in modo significativo:
  # differiscono solo tra loro, il WT sta in mezzo
  p[!s_ko & s_eeko & !s_eewt]                                     <- "KO vs EE only (WT in between)"
  # solo l'EE differisce dal WT, il KO sta in mezzo
  p[!s_ko & !s_eeko & s_eewt]                                     <- "EE vs WT only (KO in between)"
  # LRT significativo ma nessun confronto a coppie lo e'
  p[!s_ko & !s_eeko & !s_eewt]                                    <- "LRT only"
  p
}

pattern_levels <- c("KO altered, rescued by EE", "KO altered, partially rescued",
                    "KO altered, EE intermediate", "KO altered, not rescued",
                    "EE-specific", "KO vs EE only (WT in between)",
                    "EE vs WT only (KO in between)", "LRT only", "Other")

# ============================================
# Analisi di un tessuto
# ============================================
run_three_groups <- function(tissue) {
  out_dir <- file.path(deseq_dir, paste0("three_groups_", tissue))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\n========== Three groups:", tissue, "==========\n")

  cd <- samples[samples$tissue == tissue, ]
  cd$group <- factor(cd$group, levels = group_levels)
  txi <- subset_txi(txi_all, cd$sample)
  
## FROM HERE IT IS GENERIC
  dds <- DESeqDataSetFromTximport(txi, colData = cd, design = ~ group)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  # I geni del cromosoma Y riflettono solo il sesso degli animali: li escludo
  dds <- dds[!chrom_of(rownames(dds)) %in% "Y", ]

  # Test LRT: il gene cambia tra almeno due gruppi?
  dds <- DESeq(dds, test = "LRT", reduced = ~ 1)
  lrt <- results(dds)
  # Contrasti a coppie dallo stesso modello (test di Wald, stesse dispersioni)
  dds_w <- nbinomWaldTest(dds)
  ko    <- results(dds_w, contrast = c("group", "Bex3KO", "WT_Control"),    alpha = padj_cut)
  ee_wt <- results(dds_w, contrast = c("group", "Bex3KO_EE", "WT_Control"), alpha = padj_cut)
  ee_ko <- results(dds_w, contrast = c("group", "Bex3KO_EE", "Bex3KO"),     alpha = padj_cut)

  # ---------- Tabella completa ----------
  tab <- data.frame(
    gene_id   = rownames(lrt),
    gene_name = symbol_of(rownames(lrt)),
    chrom     = chrom_of(rownames(lrt)),
    baseMean  = lrt$baseMean,
    LRT_padj  = lrt$padj,
    # Per ogni confronto: log2FC, errore standard, statistica di Wald (usata per la GSEA), p-value, padj
    log2FC_Bex3KO_vs_WT        = ko$log2FoldChange,    lfcSE_Bex3KO_vs_WT        = ko$lfcSE,
    stat_Bex3KO_vs_WT          = ko$stat,              pvalue_Bex3KO_vs_WT       = ko$pvalue,
    padj_Bex3KO_vs_WT          = ko$padj,
    log2FC_Bex3KO_EE_vs_WT     = ee_wt$log2FoldChange, lfcSE_Bex3KO_EE_vs_WT     = ee_wt$lfcSE,
    stat_Bex3KO_EE_vs_WT       = ee_wt$stat,           pvalue_Bex3KO_EE_vs_WT    = ee_wt$pvalue,
    padj_Bex3KO_EE_vs_WT       = ee_wt$padj,
    log2FC_Bex3KO_EE_vs_Bex3KO = ee_ko$log2FoldChange, lfcSE_Bex3KO_EE_vs_Bex3KO = ee_ko$lfcSE,
    stat_Bex3KO_EE_vs_Bex3KO   = ee_ko$stat,           pvalue_Bex3KO_EE_vs_Bex3KO = ee_ko$pvalue,
    padj_Bex3KO_EE_vs_Bex3KO   = ee_ko$padj
  )
  tab$pattern <- ifelse(!is.na(tab$LRT_padj) & tab$LRT_padj < padj_cut,
                        classify_pattern(ko, ee_wt, ee_ko), "Not significant (LRT)")
  tab <- tab[order(tab$LRT_padj), ]
  write.csv(tab, file.path(out_dir, "three_groups_results.csv"), row.names = FALSE)

  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$gene_name <- symbol_of(rownames(norm_counts))
  write.csv(norm_counts, file.path(out_dir, "normalized_counts.csv"))

  lrt_sig <- tab[!is.na(tab$LRT_padj) & tab$LRT_padj < padj_cut, ]
  pattern_counts <- table(factor(lrt_sig$pattern, levels = pattern_levels))
  cat("Genes changing among the 3 groups (LRT padj <", padj_cut, "):", nrow(lrt_sig), "\n")
  print(pattern_counts)
  capture.output(
    cat("Genes changing among the 3 groups (LRT padj <", padj_cut, "):", nrow(lrt_sig), "\n\n"),
    print(pattern_counts),
    cat("\nPairwise contrasts (padj <", padj_cut, "):\n"),
    cat("  Bex3KO vs WT:        ", sum(ko$padj < padj_cut, na.rm = TRUE), "\n"),
    cat("  Bex3KO EE vs WT:     ", sum(ee_wt$padj < padj_cut, na.rm = TRUE), "\n"),
    cat("  Bex3KO EE vs Bex3KO: ", sum(ee_ko$padj < padj_cut, na.rm = TRUE), "\n"),
    file = file.path(out_dir, "summary.txt")
  )
  # Una lista di geni per pattern (ID Ensembl e simboli)
  for (p in setdiff(pattern_levels, "Other")) {
    sub <- lrt_sig[lrt_sig$pattern == p, ]
    fname <- gsub("[ ,]+", "_", tolower(p))
    writeLines(sub$gene_id, file.path(out_dir, paste0("genes_", fname, ".txt")))
    writeLines(unique(sub$gene_name), file.path(out_dir, paste0("genes_", fname, "_symbols.txt")))
  }
  writeLines(tab$gene_id[!is.na(tab$LRT_padj)], file.path(out_dir, "genes_background.txt"))

  vsd <- vst(dds, blind = TRUE)

  # ---------- PCA ----------
  pca <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
  pct <- round(100 * attr(pca, "percentVar"))
  pca$label <- gsub("_", " ", pca$name)
  p_pca <- ggplot(pca, aes(PC1, PC2, color = group, label = label)) +
    geom_point(size = 4) +
    geom_text(size = 3, vjust = -1.2, show.legend = FALSE) +
    scale_color_manual(values = group_colors, labels = gsub("_", " ", group_levels)) +
    xlab(paste0("PC1: ", pct[1], "% variance")) +
    ylab(paste0("PC2: ", pct[2], "% variance")) +
    labs(title = paste0("PCA: WT Control vs Bex3KO vs Bex3KO EE (", tissue, ")"), color = "Group") +
    theme_bw()
  ggsave(file.path(out_dir, "pca.png"), p_pca, width = 7.5, height = 5)

  # Il diagramma di Venn e il barplot dei DEG si creano con three_groups_venn.R e
  # three_groups_barplot.R, che leggono three_groups_results.csv

  # ---------- Heatmap dei geni LRT, ordinati per pattern ----------
  hm_genes <- lrt_sig[!lrt_sig$pattern %in% c("Other", "LRT only"), ]
  if (nrow(hm_genes) > 0) {
    mat <- assay(vsd)[hm_genes$gene_id, , drop = FALSE]
    z <- t(scale(t(mat)))
    ord_samples <- cd$sample[order(cd$group)]
    ord_genes <- hm_genes$gene_id[order(factor(hm_genes$pattern, levels = pattern_levels),
                                        -hm_genes$log2FC_Bex3KO_vs_WT)]
    hm <- data.frame(
      gene    = rep(ord_genes, times = length(ord_samples)),
      sample  = rep(ord_samples, each = length(ord_genes)),
      z       = as.vector(z[ord_genes, ord_samples]),
      pattern = rep(hm_genes$pattern[match(ord_genes, hm_genes$gene_id)], times = length(ord_samples))
    )
    hm$gene    <- factor(hm$gene, levels = rev(ord_genes))
    hm$sample  <- factor(gsub("_", " ", hm$sample), levels = gsub("_", " ", ord_samples))
    hm$pattern <- factor(hm$pattern, levels = pattern_levels)
    hm$z <- pmax(pmin(hm$z, 2), -2)
    p_hm <- ggplot(hm, aes(sample, gene, fill = z)) +
      geom_tile() +
      facet_grid(pattern ~ ., scales = "free_y", space = "free_y",
                 labeller = label_wrap_gen(width = 16)) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", name = "z-score") +
      labs(title = paste0("Genes changing among groups, by pattern (", tissue, ")"), x = NULL, y = NULL) +
      theme_minimal() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text.y = element_text(angle = 0, hjust = 0), panel.grid = element_blank())
    ggsave(file.path(out_dir, "heatmap_patterns.png"), p_hm, width = 8, height = 9, bg = "white")

    # ---------- Andamento medio di ogni pattern (z-score medio per gruppo) ----------
    zm <- sapply(group_levels, function(g) rowMeans(z[, cd$sample[cd$group == g], drop = FALSE]))
    traj <- data.frame(
      gene    = rep(rownames(zm), times = 3),
      group   = factor(rep(group_levels, each = nrow(zm)), levels = group_levels),
      z       = as.vector(zm),
      pattern = factor(rep(hm_genes$pattern[match(rownames(zm), hm_genes$gene_id)], times = 3),
                       levels = pattern_levels)
    )
    # Up/Down: per i pattern con effetto KO uso KO vs WT, altrimenti la direzione dell'EE
    hg <- hm_genes[match(rownames(zm), hm_genes$gene_id), ]
    ref_lfc <- ifelse(grepl("^KO altered", hg$pattern), hg$log2FC_Bex3KO_vs_WT,
               ifelse(hg$pattern == "KO vs EE only (WT in between)", hg$log2FC_Bex3KO_EE_vs_Bex3KO,
                      hg$log2FC_Bex3KO_EE_vs_WT))
    traj$direction <- rep(ifelse(ref_lfc > 0, "Up", "Down"), times = 3)
    n_lab <- aggregate(gene ~ pattern + direction, data = traj[traj$group == "WT_Control", ], FUN = length)
    n_lab$label <- paste0("n = ", n_lab$gene)
    p_traj <- ggplot(traj, aes(group, z, group = gene)) +
      geom_line(alpha = 0.08, color = "grey40") +
      stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 1.3, color = "black") +
      stat_summary(aes(group = 1, color = group), fun = mean, geom = "point", size = 3) +
      geom_text(data = n_lab, aes(x = 2, y = 2.3, label = label), inherit.aes = FALSE, size = 3) +
      facet_grid(direction ~ pattern, labeller = labeller(pattern = label_wrap_gen(width = 16))) +
      scale_color_manual(values = group_colors, guide = "none") +
      scale_x_discrete(labels = c("WT", "KO", "KO EE")) +
      coord_cartesian(ylim = c(-2.5, 2.5)) +
      labs(title = paste0("Expression patterns across groups (", tissue, ")"),
           x = NULL, y = "Mean z-score (VST)") +
      theme_bw()
    ggsave(file.path(out_dir, "pattern_profiles.png"), p_traj, width = 14, height = 5.5)
  }

  # ---------- Boxplot dei geni piu' significativi per i pattern principali ----------
  top <- do.call(rbind, lapply(c("KO altered, rescued by EE", "KO altered, not rescued", "EE-specific",
                                 "KO vs EE only (WT in between)"),
                               function(p) head(lrt_sig[lrt_sig$pattern == p, ], 4)))
  # (pattern principali: rescued, not rescued, EE-specific)
  if (!is.null(top) && nrow(top) > 0) {
    nc <- counts(dds, normalized = TRUE)[top$gene_id, , drop = FALSE]
    bx <- data.frame(
      gene    = rep(top$gene_name, times = ncol(nc)),
      pattern = rep(top$pattern, times = ncol(nc)),
      sample  = rep(colnames(nc), each = nrow(nc)),
      count   = as.vector(nc)
    )
    bx$group <- factor(cd$group[match(bx$sample, cd$sample)], levels = group_levels)
    bx$panel <- factor(paste0(bx$gene, "\n", bx$pattern), levels = unique(paste0(top$gene_name, "\n", top$pattern)))
    p_bx <- ggplot(bx, aes(group, count, color = group)) +
      geom_boxplot(outlier.shape = NA, color = "grey60") +
      geom_point(size = 2.5) +
      facet_wrap(~ panel, scales = "free_y", ncol = 4) +
      scale_color_manual(values = group_colors, guide = "none") +
      scale_x_discrete(labels = c("WT", "KO", "KO EE")) +
      labs(title = paste0("Top genes per pattern (", tissue, ")"), x = NULL, y = "Normalized counts") +
      theme_bw() + theme(strip.text = element_text(size = 8))
    ggsave(file.path(out_dir, "top_genes_boxplots.png"), p_bx, width = 11, height = 10)
  }

  cat("Results and plots in:", out_dir, "\n")
}

# Prima il cervelletto, poi l'ippocampo
for (tissue in c("Cb", "Hc")) run_three_groups(tissue)
