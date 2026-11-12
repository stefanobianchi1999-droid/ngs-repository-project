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
gtf_file   <- "data/reference/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
deseq_dir  <- "results/rna-seq/Exercise_GF/deseq2"

# ============================================
# Step 1: tabella dei campioni (dal modulo HaploX "Sample Information Collection Form")
#   NP1-3   Bex3KO_Ct_Hc    NP4-6   Bex3KO_Ct_Cb    NP7-9   Bex3KO_EE_Hc
#   NP10-12 Bex3KO_EE_Cb    NP13-15 Bex_WT_Ct_Hc    NP16-18 Bex_WT_Ct_Cb
# ============================================
samples <- data.frame(
  np     = paste0("NP", 1:18),
  group  = rep(c("Bex3KO", "Bex3KO", "Bex3KO_EE", "Bex3KO_EE", "WT_Control", "WT_Control"), each = 3),
  tissue = rep(c("Hc", "Cb", "Hc", "Cb", "Hc", "Cb"), each = 3),
  rep    = rep(1:3, times = 6)
)
samples$condition <- paste(samples$group, samples$tissue, sep = "_")
# Nome richiesto dal professore (es. WT_Control_Hc_1): usato in tabelle e grafici
samples$sample <- paste(samples$condition, samples$rep, sep = "_")
rownames(samples) <- samples$sample

# ============================================
# Step 2: tabella trascritto -> gene (tx2gene) dal GTF
# ============================================
# quote = "" tiene le virgolette degli attributi (gene_id "ENSMUSG...";)
gtf <- read.delim(gzfile(gtf_file), header = FALSE, comment.char = "#", quote = "",
                  colClasses = c("character", "NULL", "character", rep("NULL", 5), "character"))
colnames(gtf) <- c("chrom", "feature", "attributes")
gtf <- gtf[gtf$feature == "transcript", ]

# Estrae il valore di un attributo GTF (es. gene_id "ENSMUSG..."; -> ENSMUSG...)
get_attr <- function(attrs, key) {
  sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", attrs)
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
# Step 3: import di tutti i quant.sf con tximport (somma trascritti -> geni)
# ============================================
files <- file.path(salmon_dir, samples$np, "quant.sf")
names(files) <- samples$sample
# quant.sf ha gli ID con versione (ENSMUST...1), il GTF senza -> ignoreTxVersion
txi_all <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Sottoinsieme di txi con solo alcuni campioni
subset_txi <- function(txi, keep) {
  txi$abundance <- txi$abundance[, keep, drop = FALSE]
  txi$counts    <- txi$counts[, keep, drop = FALSE]
  txi$length    <- txi$length[, keep, drop = FALSE]
  txi
}

