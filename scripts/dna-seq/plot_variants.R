library(ggplot2)

# Read VCF, skipping header lines (##...)
vcf_lines <- readLines("results/dna-seq/variants/SRR2584866.vcf.gz")
vcf_lines <- vcf_lines[!grepl("^##", vcf_lines)]

# Parse into a data frame
vcf <- read.table(text = vcf_lines, header = TRUE, comment.char = "", sep = "\t")
colnames(vcf)[1] <- "CHROM"

# Classify variant type: SNP vs INDEL
vcf$type <- ifelse(nchar(vcf$REF) == 1 & nchar(vcf$ALT) == 1, "SNP", "INDEL")

# Plot 1: distribution of variants along the genome
p1 <- ggplot(vcf, aes(x = POS, fill = type)) +
  geom_histogram(binwidth = 50000, color = "white") +
  scale_fill_manual(values = c("SNP" = "steelblue", "INDEL" = "firebrick")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Distribution of Variants Along the E. coli Genome",
    subtitle = paste0("SRR2584866 (generation 50,000) vs REL606 reference — ", nrow(vcf), " variants total"),
    x = "Genomic position (bp)",
    y = "Number of variants per 50kb bin",
    fill = "Variant type"
  )

ggsave("results/dna-seq/variant_distribution.png", plot = p1, width = 10, height = 6, dpi = 300)

# Plot 2: SNP vs INDEL counts
p2 <- ggplot(vcf, aes(x = type, fill = type)) +
  geom_bar() +
  scale_fill_manual(values = c("SNP" = "steelblue", "INDEL" = "firebrick")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Variant Types",
    x = "", y = "Count"
  ) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5)

ggsave("results/dna-seq/variant_types.png", plot = p2, width = 6, height = 6, dpi = 300)

cat("Variant type breakdown:\n")
print(table(vcf$type))

