# Load DESeq2
library(DESeq2)

# Read the featureCounts output, skipping the comment line
counts_raw <- read.table("results/counts/gene_counts.txt", header = TRUE, skip = 1, row.names = 1)

# Keep only the count columns (drop Chr, Start, End, Strand, Length)
counts <- counts_raw[, 6:ncol(counts_raw)]

# Clean up column names (remove path and file suffix)
colnames(counts) <- gsub("results.aligned\\.", "", colnames(counts))
colnames(counts) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(counts))

# Define sample conditions (must match column order!)
condition <- factor(c("normal", "normal", "normal", "tumor", "tumor", "tumor"))
coldata <- data.frame(row.names = colnames(counts), condition = condition)

# Build DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Order by adjusted p-value (most significant genes first)
res <- res[order(res$padj), ]

# Save results to CSV
write.csv(as.data.frame(res), file = "results/counts/deseq2_results.csv")

# Print summary
summary(res)

cat("\nTop 10 differentially expressed genes:\n")
print(head(res, 10))
