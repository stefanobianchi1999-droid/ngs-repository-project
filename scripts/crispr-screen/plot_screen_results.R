library(ggplot2)
library(dplyr)

# ============================================================
# Load gene-level summary (from mageck test)
# ============================================================
gene_res <- read.delim("results/crispr-screen/counts/screen_test.gene_summary.txt", sep = "\t")

# Use the "neg" direction (genes depleted after selection = important for survival)
gene_res$fdr_use <- gene_res$neg.fdr
gene_res$lfc_use <- gene_res$neg.lfc
gene_res$pval_use <- gene_res$neg.p.value

gene_res$significant <- ifelse(gene_res$fdr_use < 0.05, "Significant (FDR<0.05)", "Not significant")

# ============================================================
# 1. Volcano plot
# ============================================================
volcano <- ggplot(gene_res, aes(x = lfc_use, y = -log10(pval_use), color = significant)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Significant (FDR<0.05)" = "firebrick", "Not significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    title = "CRISPRi Screen: Gene-level Volcano Plot",
    subtitle = "E. coli tiling library, MOPS selection vs initial",
    x = "log2 Fold Change (guide abundance)",
    y = "-log10(p-value)",
    color = "Category"
  )
ggsave("results/crispr-screen/volcano_plot.png", plot = volcano, width = 8, height = 6, dpi = 300)

# ============================================================
# 2. Rank plot (waterfall)
# ============================================================
gene_res_sorted <- gene_res[order(gene_res$lfc_use), ]
gene_res_sorted$rank <- seq_len(nrow(gene_res_sorted))

rankplot <- ggplot(gene_res_sorted, aes(x = rank, y = lfc_use, color = significant)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Significant (FDR<0.05)" = "firebrick", "Not significant" = "grey70")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  labs(
    title = "CRISPRi Screen: Gene Rank Plot",
    subtitle = "Genes ranked by fold change (most depleted on the left)",
    x = "Gene rank",
    y = "log2 Fold Change",
    color = "Category"
  )
ggsave("results/crispr-screen/rank_plot.png", plot = rankplot, width = 8, height = 6, dpi = 300)

# ============================================================
# 3. Top hits bar chart
# ============================================================
top_hits <- gene_res %>% arrange(pval_use) %>% head(15)
top_hits$id <- factor(top_hits$id, levels = rev(top_hits$id))

barplot_top <- ggplot(top_hits, aes(x = id, y = lfc_use, fill = lfc_use < 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick"), guide = "none") +
  theme_bw(base_size = 14) +
  labs(
    title = "Top 15 Hit Genes",
    subtitle = "E. coli CRISPRi screen, MOPS vs initial",
    x = "Gene",
    y = "log2 Fold Change"
  )
ggsave("results/crispr-screen/top_hits_barplot.png", plot = barplot_top, width = 8, height = 6, dpi = 300)

# ============================================================
# 4. sgRNA-level scatter plot (initial vs mops abundance)
# ============================================================
sgrna_counts <- read.delim("results/crispr-screen/counts/screen.count.txt", sep = "\t")

# log-transform counts (add 1 to avoid log(0))
sgrna_counts$log_initial <- log2(sgrna_counts$initial + 1)
sgrna_counts$log_mops <- log2(sgrna_counts$mops + 1)

scatter <- ggplot(sgrna_counts, aes(x = log_initial, y = log_mops)) +
  geom_point(alpha = 0.4, size = 1.5, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 14) +
  labs(
    title = "sgRNA Abundance: Initial vs MOPS",
    subtitle = "Each point is one sgRNA; dashed line = no change",
    x = "log2(read count) - Initial",
    y = "log2(read count) - MOPS"
  )
ggsave("results/crispr-screen/sgrna_scatter.png", plot = scatter, width = 8, height = 6, dpi = 300)

cat("All 4 plots saved to results/crispr-screen/\n")
