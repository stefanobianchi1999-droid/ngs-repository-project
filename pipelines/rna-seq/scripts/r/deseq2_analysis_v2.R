source("pipelines/rna-seq/scripts/r/prep_salmon.R")
source("pipelines/rna-seq/scripts/r/deseq2_charts.R")

# tissues <- snakemake@params[["tissues"]]
#
# prep <- prep_salmon(
#   groups     = snakemake@params[["groups"]],
#   reference  = snakemake@params[["reference"]],
#   tissues    = tissues,
#   reps       = snakemake@params[["reps"]],
#   gtf_file   = snakemake@params[["gtf_file"]],
#   sample_dir = snakemake@params[["sample_dir"]]
# )

# --- Test senza Snakemake: valori costanti presi da deseq2_pipeline/config.yml ---
tissues <- c("Hc", "Cb")

prep <- prep_salmon(
  groups     = c("Bex3KO", "Bex3KO_EE", "WT_Control"),   # ordine delle cartelle NP
  reference  = "WT_Control",
  tissues    = tissues,
  reps       = 3,
  gtf_file   = "data/reference/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz",
  sample_dir = "results/rna-seq/Exercise_GF/salmon"
)

# Unpacks prep into plain variables so the generic code below can keep
# using their names.
txi_all <- prep$txi_all
samples <- prep$samples
subset_txi <- prep$subset_txi
symbol_of <- prep$symbol_of
chrom_of <- prep$chrom_of
# Ordine di analisi (controllo per primo), NON l'ordine delle cartelle NP:
# lo usano combn per la direzione dei contrasti e i grafici
groups <- prep$group_levels
group_colors <- setNames(hcl.colors(length(groups), "Dark 3"), groups)

#rm(prep)


# ============================================
# Classificazione dei pattern a partire dai tre contrasti a coppie
#   ko = Bex3KO vs WT, ee_wt = Bex3KO EE vs WT, ee_ko = Bex3KO EE vs Bex3KO
# ============================================

# padj_cut     <- snakemake@params[["padj_cut"]]
padj_cut     <- 0.05

classify_pattern <- function(ko, ee_wt, ee_ko) {
  sig <- function(x) !is.na(x$padj) & x$padj < padj_cut
  s_ko <- sig(ko)
  s_eewt <- sig(ee_wt)
  s_eeko <- sig(ee_ko)

  d_ko <- sign(ko$log2FoldChange)
  d_eewt <- sign(ee_wt$log2FoldChange)
  d_eeko <- sign(ee_ko$log2FoldChange)

  p <- rep("Other", length(s_ko))
  # alterato nel KO, l'EE lo riporta verso il WT
  p[s_ko & s_eeko & d_eeko == -d_ko & !s_eewt] <- "KO altered, rescued by EE"
  p[s_ko & s_eeko & d_eeko == -d_ko & s_eewt & d_eewt == d_ko] <- "KO altered, partially rescued"
  # alterato nel KO, l'EE non lo corregge
  p[s_ko & !s_eeko & s_eewt & d_eewt == d_ko] <- "KO altered, not rescued"
  # alterato nel KO, l'EE sta a meta' (nessun contrasto con l'EE significativo)
  p[s_ko & !s_eeko & !s_eewt] <- "KO altered, EE intermediate"
  # non alterato nel KO, cambiato solo dall'EE
  p[!s_ko & s_eeko & s_eewt & d_eeko == d_eewt] <- "EE-specific"
  # KO ed EE si discostano dal WT in direzioni opposte, nessuno dei due in modo significativo:
  # differiscono solo tra loro, il WT sta in mezzo
  p[!s_ko & s_eeko & !s_eewt] <- "KO vs EE only (WT in between)"
  # solo l'EE differisce dal WT, il KO sta in mezzo
  p[!s_ko & !s_eeko & s_eewt] <- "EE vs WT only (KO in between)"
  # LRT significativo ma nessun confronto a coppie lo e'
  p[!s_ko & !s_eeko & !s_eewt] <- "LRT only"
  p
}

pattern_levels <- c("KO altered, rescued by EE", "KO altered, partially rescued",
                    "KO altered, EE intermediate", "KO altered, not rescued",
                    "EE-specific", "KO vs EE only (WT in between)",
                    "EE vs WT only (KO in between)", "LRT only", "Other")

# ============================================
# Analisi di un tessuto
# ============================================
# output_dir  <- snakemake@params[["output_dir"]]
output_dir  <- "results/rna-seq/Exercise_GF/deseq2"

run_three_groups <- function(prep, tissue) {
  subfolder_name <- "three_groups_"
  out_dir <- file.path(output_dir, paste0(subfolder_name, tissue))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\n========== Three groups:", tissue, "==========\n")

  dds <- dds_salmon_creation(prep, tissue)
  cd  <- as.data.frame(colData(dds))
  lrt <- results(dds)

  # Contrasti a coppie dallo stesso modello (test di Wald, stesse dispersioni)
  dds_w <- nbinomWaldTest(dds)

  # Tutte le coppie di gruppi: c(numeratore, denominatore), il denominatore e' il gruppo che viene prima
  # in groups (= group_levels: reference, poi gli altri nell'ordine delle cartelle NP)
  comparisons <- combn(groups, 2, FUN = function(p) c(p[2], p[1]), simplify = FALSE)
  names(comparisons) <- vapply(comparisons, paste, character(1), collapse = "_vs_")

  res <- lapply(comparisons, function(cmp) {
    results(dds_w, contrast = c("group", cmp[1], cmp[2]), alpha = padj_cut)
  })

  block <- Map(function(result, nm) {
    df <- data.frame(
      log2FC = result$log2FoldChange,
      lfcSE  = result$lfcSE,
      stat   = result$stat,
      pvalue = result$pvalue,
      padj   = result$padj
    )
    names(df) <- paste0(names(df), "_", nm)   # es. log2FC_Bex3KO_vs_WT_Control
    df
  }, res, names(res))
  block_df <- do.call(cbind, unname(block))

  tab <- data.frame(
    gene_id   = rownames(lrt),
    gene_name = symbol_of(rownames(lrt)),
    chrom     = chrom_of(rownames(lrt)),
    baseMean  = lrt$baseMean,
    LRT_padj  = lrt$padj
  )
  tab <- cbind(tab, block_df)


  #tanti saluti li abbiamo convertiti
  ko    <- results(dds_w, contrast = c("group", "Bex3KO", "WT_Control"),    alpha = padj_cut)
  ee_wt <- results(dds_w, contrast = c("group", "Bex3KO_EE", "WT_Control"), alpha = padj_cut)
  ee_ko <- results(dds_w, contrast = c("group", "Bex3KO_EE", "Bex3KO"),     alpha = padj_cut)
  
  tab$pattern <- ifelse(!is.na(tab$LRT_padj) & tab$LRT_padj < padj_cut,
                        classify_pattern(ko, ee_wt, ee_ko), "Not significant (LRT)")
  
  tab <- tab[order(tab$LRT_padj), ]

  file_path_three_groups <- file.path(out_dir, "three_groups_results.csv")
  write.csv(tab, file_path_three_groups, row.names = FALSE)

  file_path_normalized <- file.path(out_dir, "normalized_counts.csv")
  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$gene_name <- symbol_of(rownames(norm_counts))
  write.csv(norm_counts, file_path_normalized)

  lrt_sig <- tab[!is.na(tab$LRT_padj) & tab$LRT_padj < padj_cut, ]
  pattern_counts <- table(factor(lrt_sig$pattern, levels = pattern_levels))

  output_texts <- Map(function(result, nm) {
    paste(nm, sum(result$padj < padj_cut, na.rm = TRUE))
  }, res, names(res))
  
  file_path_summary <- file.path(out_dir, "summary.txt")

  capture.output(
    cat("Genes changing among the", length(groups), "groups (LRT padj <", padj_cut, "):", nrow(lrt_sig), "\n\n"),
    print(pattern_counts),
    cat("\nPairwise contrasts (padj <", padj_cut, "):\n"),
    cat(unlist(output_texts), sep = "\n"),
    file = file_path_summary
  )

  file_path_gene_background <- file.path(out_dir, "genes_background.txt")
  # Una lista di geni per pattern (ID Ensembl e simboli)
  for (p in setdiff(pattern_levels, "Other")) {
    sub <- lrt_sig[lrt_sig$pattern == p, ]
    fname <- gsub("[ ,]+", "_", tolower(p))
    writeLines(sub$gene_id,file.path(out_dir, paste0("genes_", fname, ".txt")))
    writeLines(unique(sub$gene_name), file.path(out_dir, paste0("genes_", fname, "_symbols.txt")))
  }
  writeLines(tab$gene_id[!is.na(tab$LRT_padj)], file_path_gene_background)

  vsd <- vst(dds, blind = TRUE)

  # ---------- PCA ----------
  output_file <- file.path(out_dir, "pca.png")
  PCA_chart_drawer(vsd = vsd, groups = groups, tissue = tissue,
                   group_colors = group_colors, output_file = output_file)

  # # Il diagramma di Venn e il barplot dei DEG si creano con three_groups_venn.R e
  # # three_groups_barplot.R, che leggono three_groups_results.csv

  # # ---------- Heatmap dei geni LRT, ordinati per pattern ----------
  # hm_genes <- lrt_sig[!lrt_sig$pattern %in% c("Other", "LRT only"), ]
  # if (nrow(hm_genes) > 0) {
  #   mat <- assay(vsd)[hm_genes$gene_id, , drop = FALSE]
  #   z <- t(scale(t(mat)))
  #   ord_samples <- cd$sample[order(cd$group)]
  #   ord_genes <- hm_genes$gene_id[order(factor(hm_genes$pattern, levels = pattern_levels),
  #                                       -hm_genes$log2FC_Bex3KO_vs_WT)]
  #   hm <- data.frame(
  #     gene    = rep(ord_genes, times = length(ord_samples)),
  #     sample  = rep(ord_samples, each = length(ord_genes)),
  #     z       = as.vector(z[ord_genes, ord_samples]),
  #     pattern = rep(hm_genes$pattern[match(ord_genes, hm_genes$gene_id)], times = length(ord_samples))
  #   )
  #   hm$gene    <- factor(hm$gene, levels = rev(ord_genes))
  #   hm$sample  <- factor(gsub("_", " ", hm$sample), levels = gsub("_", " ", ord_samples))
  #   hm$pattern <- factor(hm$pattern, levels = pattern_levels)
  #   hm$z <- pmax(pmin(hm$z, 2), -2)
  #   p_hm <- ggplot(hm, aes(sample, gene, fill = z)) +
  #     geom_tile() +
  #     facet_grid(pattern ~ ., scales = "free_y", space = "free_y",
  #                labeller = label_wrap_gen(width = 16)) +
  #     scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", name = "z-score") +
  #     labs(title = paste0("Genes changing among groups, by pattern (", tissue, ")"), x = NULL, y = NULL) +
  #     theme_minimal() +
  #     theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
  #           axis.text.x = element_text(angle = 45, hjust = 1),
  #           strip.text.y = element_text(angle = 0, hjust = 0), panel.grid = element_blank())
  #   ggsave(file.path(out_dir, "heatmap_patterns.png"), p_hm, width = 8, height = 9, bg = "white")

  #   # ---------- Andamento medio di ogni pattern (z-score medio per gruppo) ----------
  #   zm <- sapply(groups, function(g) rowMeans(z[, cd$sample[cd$group == g], drop = FALSE]))
  #   traj <- data.frame(
  #     gene    = rep(rownames(zm), times = 3),
  #     group   = factor(rep(groups, each = nrow(zm)), levels = groups),
  #     z       = as.vector(zm),
  #     pattern = factor(rep(hm_genes$pattern[match(rownames(zm), hm_genes$gene_id)], times = 3),
  #                      levels = pattern_levels)
  #   )
  #   # Up/Down: per i pattern con effetto KO uso KO vs WT, altrimenti la direzione dell'EE
  #   hg <- hm_genes[match(rownames(zm), hm_genes$gene_id), ]
  #   ref_lfc <- ifelse(grepl("^KO altered", hg$pattern), hg$log2FC_Bex3KO_vs_WT,
  #              ifelse(hg$pattern == "KO vs EE only (WT in between)", hg$log2FC_Bex3KO_EE_vs_Bex3KO,
  #                     hg$log2FC_Bex3KO_EE_vs_WT))
  #   traj$direction <- rep(ifelse(ref_lfc > 0, "Up", "Down"), times = 3)
  #   n_lab <- aggregate(gene ~ pattern + direction, data = traj[traj$group == "WT_Control", ], FUN = length)
  #   n_lab$label <- paste0("n = ", n_lab$gene)
  #   p_traj <- ggplot(traj, aes(group, z, group = gene)) +
  #     geom_line(alpha = 0.08, color = "grey40") +
  #     stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 1.3, color = "black") +
  #     stat_summary(aes(group = 1, color = group), fun = mean, geom = "point", size = 3) +
  #     geom_text(data = n_lab, aes(x = 2, y = 2.3, label = label), inherit.aes = FALSE, size = 3) +
  #     facet_grid(direction ~ pattern, labeller = labeller(pattern = label_wrap_gen(width = 16))) +
  #     scale_color_manual(values = group_colors, guide = "none") +
  #     scale_x_discrete(labels = c("WT", "KO", "KO EE")) +
  #     coord_cartesian(ylim = c(-2.5, 2.5)) +
  #     labs(title = paste0("Expression patterns across groups (", tissue, ")"),
  #          x = NULL, y = "Mean z-score (VST)") +
  #     theme_bw()
  #   ggsave(file.path(out_dir, "pattern_profiles.png"), p_traj, width = 14, height = 5.5)
  # }

  # # ---------- Boxplot dei geni piu' significativi per i pattern principali ----------
  # top <- do.call(rbind, lapply(c("KO altered, rescued by EE", "KO altered, not rescued", "EE-specific",
  #                                "KO vs EE only (WT in between)"),
  #                              function(p) head(lrt_sig[lrt_sig$pattern == p, ], 4)))
  # # (pattern principali: rescued, not rescued, EE-specific)
  # if (!is.null(top) && nrow(top) > 0) {
  #   nc <- counts(dds, normalized = TRUE)[top$gene_id, , drop = FALSE]
  #   bx <- data.frame(
  #     gene    = rep(top$gene_name, times = ncol(nc)),
  #     pattern = rep(top$pattern, times = ncol(nc)),
  #     sample  = rep(colnames(nc), each = nrow(nc)),
  #     count   = as.vector(nc)
  #   )
  #   bx$group <- factor(cd$group[match(bx$sample, cd$sample)], levels = groups)
  #   bx$panel <- factor(paste0(bx$gene, "\n", bx$pattern), levels = unique(paste0(top$gene_name, "\n", top$pattern)))
  #   p_bx <- ggplot(bx, aes(group, count, color = group)) +
  #     geom_boxplot(outlier.shape = NA, color = "grey60") +
  #     geom_point(size = 2.5) +
  #     facet_wrap(~ panel, scales = "free_y", ncol = 4) +
  #     scale_color_manual(values = group_colors, guide = "none") +
  #     scale_x_discrete(labels = c("WT", "KO", "KO EE")) +
  #     labs(title = paste0("Top genes per pattern (", tissue, ")"), x = NULL, y = "Normalized counts") +
  #     theme_bw() + theme(strip.text = element_text(size = 8))
  #   ggsave(file.path(out_dir, "top_genes_boxplots.png"), p_bx, width = 11, height = 10)
  # }

  # cat("Results and plots in:", out_dir, "\n")
}

# Prima il cervelletto, poi l'ippocampo
for (tissue in tissues) run_three_groups(prep, tissue)
