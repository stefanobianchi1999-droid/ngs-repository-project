# Load DESeq2 results
res <- read.csv("results/counts/deseq2_results.csv", row.names = 1)

# Load gene ID -> gene name mapping (built from the GTF file)
gene_names <- read.table("data/reference/gene_id_to_name.tsv",
                          sep = "\t", header = FALSE,
                          col.names = c("gene_id", "gene_name"))

# Match gene names to results by row name (Ensembl ID)
res$gene_name <- gene_names$gene_name[match(rownames(res), gene_names$gene_id)]

# Reorder columns so gene_name appears first, right after the row names
res <- res[, c("gene_name", setdiff(colnames(res), "gene_name"))]

# Keep only significant genes, sorted by adjusted p-value
res_sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
res_sig <- res_sig[order(res_sig$padj), ]

# Save top 20 significant genes with real names
top20 <- head(res_sig, 20)
write.csv(top20, file = "results/counts/top20_genes_annotated.csv")

cat("Top 20 significant genes (annotated):\n")
print(top20[, c("gene_name", "log2FoldChange", "padj")])
