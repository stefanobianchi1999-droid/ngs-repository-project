# ============================================
# Exercise_GF: Salmon -> tximport -> DESeq2 -> grafici
# Topo Bex3 KO, GRCm38 / Ensembl 102. Ogni confronto ha la sua sottocartella
# in results/rna-seq/Exercise_GF/deseq2/<confronto>/
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/salmon_deseq2_exercise_gf.R")
#
# Per aggiungere campioni: aggiorna la tabella `samples` (Step 1).
# Per aggiungere confronti: aggiungi una riga in `comparisons` (Step 4).
# ============================================
source("pipelines/rna-seq/scripts/r/exercise_gf_setup.R")

# ============================================
# Step 4: confronti da fare (test = gruppo di interesse, ref = riferimento)
# log2FoldChange > 0 = piu' espresso in `test`
# Direzioni come nel modulo HaploX (Comparison Group 6-13)
# ============================================
comparisons <- data.frame(
  name = c(
    # confronti tra gruppi, separati per tessuto
    "Hc_Bex3KO_vs_WT_Control",     "Hc_Bex3KO_EE_vs_WT_Control",  "Hc_Bex3KO_EE_vs_Bex3KO",
    "Cb_Bex3KO_vs_WT_Control",     "Cb_Bex3KO_EE_vs_WT_Control",  "Cb_Bex3KO_EE_vs_Bex3KO",
    # confronti tra tessuti, dentro lo stesso gruppo
    "WT_Control_Hc_vs_Cb",         "Bex3KO_Hc_vs_Cb",             "Bex3KO_EE_Hc_vs_Cb"),
  test = c("Bex3KO_Hc",    "Bex3KO_EE_Hc", "Bex3KO_EE_Hc",
           "Bex3KO_Cb",    "Bex3KO_EE_Cb", "Bex3KO_EE_Cb",
           "WT_Control_Hc", "Bex3KO_Hc",   "Bex3KO_EE_Hc"),
  ref  = c("WT_Control_Hc", "WT_Control_Hc", "Bex3KO_Hc",
           "WT_Control_Cb", "WT_Control_Cb", "Bex3KO_Cb",
           "WT_Control_Cb", "Bex3KO_Cb",     "Bex3KO_EE_Cb")
)

# ============================================
# Step 5: funzione che esegue un confronto e salva risultati + grafici
# ============================================
run_comparison <- function(name, test, ref, padj_cut = 0.05, lfc_cut = 1) {
  out_dir <- file.path(deseq_dir, name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\n========== ", name, ": ", test, " vs ", ref, " ==========\n", sep = "")

  # DESeq2 solo sui campioni dei due gruppi confrontati
  cd <- samples[samples$condition %in% c(test, ref), ]
  cd$condition <- factor(cd$condition, levels = c(ref, test))
  txi <- subset_txi(txi_all, cd$sample)

  dds <- DESeqDataSetFromTximport(txi, colData = cd, design = ~ condition)
  # Tolgo i geni quasi mai visti (< 10 read in almeno 3 campioni)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", test, ref), alpha = padj_cut)
  summary(res)
  capture.output(summary(res), file = file.path(out_dir, "deseq2_summary.txt"))

  res_df <- as.data.frame(res)
  res_df$gene_name <- symbol_of(rownames(res_df))
  res_df <- res_df[order(res_df$padj), c("gene_name", setdiff(colnames(res_df), "gene_name"))]
  write.csv(res_df, file.path(out_dir, "deseq2_results.csv"))

  norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_counts$gene_name <- symbol_of(rownames(norm_counts))
  write.csv(norm_counts, file.path(out_dir, "normalized_counts.csv"))

  # Liste di geni significativi, pronte per g:Profiler (ID) ed Enrichr (simboli)
  sig  <- res_df[!is.na(res_df$padj) & res_df$padj < padj_cut, ]
  up   <- sig[sig$log2FoldChange >  lfc_cut, ]
  down <- sig[sig$log2FoldChange < -lfc_cut, ]
  writeLines(rownames(up),   file.path(out_dir, paste0("genes_up_in_", test, ".txt")))
  writeLines(rownames(down), file.path(out_dir, paste0("genes_up_in_", ref, ".txt")))
  writeLines(unique(up$gene_name),   file.path(out_dir, paste0("genes_up_in_", test, "_symbols.txt")))
  writeLines(unique(down$gene_name), file.path(out_dir, paste0("genes_up_in_", ref, "_symbols.txt")))
  writeLines(rownames(res_df)[!is.na(res_df$padj)], file.path(out_dir, "genes_background.txt"))
  cat("Up in", test, ":", nrow(up), "- up in", ref, ":", nrow(down),
      "(padj <", padj_cut, ", |log2FC| >", lfc_cut, ")\n")

  # ---------- Grafici ----------
  vsd <- vst(dds, blind = TRUE)

  # PCA
  pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  pct <- round(100 * attr(pca, "percentVar"))
  pca$label <- gsub("_", " ", pca$name)
  p_pca <- ggplot(pca, aes(PC1, PC2, color = condition, label = label)) +
    geom_point(size = 4) +
    geom_text(vjust = -1, show.legend = FALSE) +
    xlab(paste0("PC1: ", pct[1], "% variance")) +
    ylab(paste0("PC2: ", pct[2], "% variance")) +
    labs(title = paste("PCA:", test, "vs", ref), color = "Condition") +
    theme_bw()
  ggsave(file.path(out_dir, "pca.png"), p_pca, width = 7, height = 5)

  # Heatmap delle distanze tra campioni
  png(file.path(out_dir, "sample_distances.png"), width = 700, height = 600)
  dist_mat <- as.matrix(dist(t(assay(vsd))))
  rownames(dist_mat) <- colnames(dist_mat) <- gsub("_", " ", colnames(dist_mat))
  heatmap(dist_mat, symm = TRUE, margins = c(12, 12), main = "Sample-to-sample distances (VST)")
  dev.off()

  # MA plot
  png(file.path(out_dir, "ma_plot.png"), width = 700, height = 600)
  plotMA(res, ylim = c(-10, 10), main = paste("MA plot:", test, "vs", ref))
  dev.off()

  # Volcano plot
  volc <- res_df[!is.na(res_df$padj), ]
  # padj = 0 (sotto la precisione del computer) darebbe -log10 = Inf: lo fisso a 1e-300
  volc$padj <- pmax(volc$padj, 1e-300)
  lab_up   <- paste("Up in", test)
  lab_down <- paste("Up in", ref)
  volc$signif <- ifelse(volc$padj < padj_cut & volc$log2FoldChange >  lfc_cut, lab_up,
                 ifelse(volc$padj < padj_cut & volc$log2FoldChange < -lfc_cut, lab_down, "Not significant"))
  top_labels <- head(volc[volc$signif != "Not significant", ], 15)
  colors <- setNames(c("firebrick", "steelblue", "grey70"), c(lab_up, lab_down, "Not significant"))
  p_volc <- ggplot(volc, aes(log2FoldChange, -log10(padj), color = signif)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_text(data = top_labels, aes(label = gene_name), size = 3, vjust = -0.7,
              show.legend = FALSE) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    labs(title = paste("Volcano plot:", test, "vs", ref),
         x = "log2 fold change", y = "-log10 adjusted p-value", color = "") +
    theme_bw()
  ggsave(file.path(out_dir, "volcano.png"), p_volc, width = 8, height = 6)

  # Heatmap dei 50 geni piu' significativi (valori VST centrati per gene)
  top50 <- head(rownames(res_df)[!is.na(res_df$padj)], 50)
  mat <- assay(vsd)[top50, ]
  mat <- mat - rowMeans(mat)
  rownames(mat) <- symbol_of(rownames(mat))
  colnames(mat) <- gsub("_", " ", colnames(mat))
  png(file.path(out_dir, "heatmap_top50.png"), width = 700, height = 1000)
  heatmap(mat, scale = "none", margins = c(12, 10), cexRow = 0.7,
          col = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
          main = "Top 50 genes by adjusted p-value")
  dev.off()

  cat("Results and plots in:", out_dir, "\n")
}

# ============================================
# Step 6: PCA di tutti i 18 campioni (panoramica generale)
# ============================================
dds_all <- DESeqDataSetFromTximport(txi_all, colData = samples, design = ~ condition)
vsd_all <- vst(dds_all, blind = TRUE)
pca_all <- plotPCA(vsd_all, intgroup = c("group", "tissue"), returnData = TRUE)
pct_all <- round(100 * attr(pca_all, "percentVar"))
# plotPCA sovrascrive la colonna "group" con gruppo:tessuto -> la ricostruisco
pca_all$Group <- gsub("_", " ", samples[pca_all$name, "group"])
pca_all$label <- gsub("_", " ", pca_all$name)
p_all <- ggplot(pca_all, aes(PC1, PC2, color = Group, shape = tissue, label = label)) +
  geom_point(size = 4) +
  geom_text(size = 2.5, vjust = -1.2, show.legend = FALSE) +
  xlab(paste0("PC1: ", pct_all[1], "% variance")) +
  ylab(paste0("PC2: ", pct_all[2], "% variance")) +
  labs(title = "PCA: all samples", color = "Group", shape = "Tissue") +
  theme_bw()
ggsave(file.path(deseq_dir, "pca_all_samples.png"), p_all, width = 7, height = 5)

# ============================================
# Step 7: eseguo tutti i confronti
# ============================================
for (i in seq_len(nrow(comparisons))) {
  run_comparison(comparisons$name[i], comparisons$test[i], comparisons$ref[i])
}
