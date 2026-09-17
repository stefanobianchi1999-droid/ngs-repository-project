library(DESeq2)
library(ggplot2)

# Load DESeq2 results (already computed and saved as CSV)
res <- read.csv("results/counts/deseq2_results.csv", row.names = 1)

# Remove genes with NA padj (DESeq2 assigns NA to very low-count genes,
# automatically filtered out and not testable)
res <- res[!is.na(res$padj), ]

# Define significance categories for coloring
res$significant <- "Not significant"
res$significant[res$padj < 0.05 & res$log2FoldChange > 1] <- "Up in tumor"
res$significant[res$padj < 0.05 & res$log2FoldChange < -1] <- "Down in tumor"

# Build the volcano plot
volcano <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Up in tumor" = "firebrick",
    "Down in tumor" = "steelblue",
    "Not significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Differential Gene Expression: Tumor vs Normal",
    subtitle = "HCC1395 cell line, chromosome 22",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Category"
  )

# Save the plot
ggsave("results/counts/volcano_plot.png", plot = volcano, width = 8, height = 6, dpi = 300)

cat("Volcano plot saved to results/counts/volcano_plot.png\n")
cat("\nGene counts by category:\n")
print(table(res$significant))
