# Load DESeq2
library(DESeq2)
library(ggplot2)

# ============================================
# Argomenti da riga di comando (passati dalla regola Snakemake `deseq2_analysis`)
#   1. counts_file   - path a gene_counts.txt (output di featureCounts)
#   2. conditions    - mapping sample->condition, formato "sample1:cond1,sample2:cond2,..."
#   3. results_csv   - path di output per la tabella dei risultati DESeq2
#   4. volcano_png   - path di output per il volcano plot
# ============================================
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
conditions_arg <- args[2]
results_csv <- args[3]
volcano_png <- args[4]

# Read the featureCounts output, skipping the comment line
counts_raw <- read.table(counts_file, header = TRUE, skip = 1, row.names = 1, check.names = FALSE)

# Keep only the count columns (drop Chr, Start, End, Strand, Length)
counts <- counts_raw[, 6:ncol(counts_raw), drop = FALSE]

# Clean up column names (remove path and file suffix) -> lasciano solo il sample id
colnames(counts) <- gsub(".*/", "", colnames(counts))
colnames(counts) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(counts))

# Parse "sample1:cond1,sample2:cond2" in una named list sample -> condition
pairs <- strsplit(strsplit(conditions_arg, ",")[[1]], ":")
condition_map <- setNames(
  vapply(pairs, `[`, character(1), 2),
  vapply(pairs, `[`, character(1), 1)
)

# Ordina le condizioni secondo l'ordine reale delle colonne nella matrice di conteggi
# (evita il rischio di un mismatch silenzioso se l'ordine dei sample cambia)
condition <- factor(condition_map[colnames(counts)])
coldata <- data.frame(row.names = colnames(counts), condition = condition)

if (nlevels(condition) < 2) {
  stop("DESeq2 richiede almeno 2 condizioni distinte (con repliche) per calcolare il contrasto. ",
       "Aggiorna 'conditions' in config.yml con piu' sample/gruppi.")
}

# Build DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Order by adjusted p-value (most significant genes first)
res <- res[order(res$padj), ]

# Save results to CSV
write.csv(as.data.frame(res), file = results_csv)

# Print summary
summary(res)

cat("\nTop 10 differentially expressed genes:\n")
print(head(res, 10))

# ============================================
# Volcano plot: log2FoldChange (asse X) vs -log10(padj) (asse Y)
# ============================================
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]

padj_threshold <- 0.05
lfc_threshold <- 1

res_df$significance <- "Not significant"
res_df$significance[res_df$padj < padj_threshold & res_df$log2FoldChange > lfc_threshold] <- "Up"
res_df$significance[res_df$padj < padj_threshold & res_df$log2FoldChange < -lfc_threshold] <- "Down"
res_df$significance <- factor(res_df$significance, levels = c("Down", "Not significant", "Up"))

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Down" = "#2166AC", "Not significant" = "grey70", "Up" = "#B2182B")) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey40") +
  labs(title = "Volcano plot", x = "log2 Fold Change", y = "-log10(padj)", color = "Significance") +
  theme_minimal()

ggsave(volcano_png, plot = volcano, width = 8, height = 6, dpi = 150)
