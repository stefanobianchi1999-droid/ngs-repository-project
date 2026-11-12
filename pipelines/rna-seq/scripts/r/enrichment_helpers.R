# ============================================
# Exercise_GF: funzioni comuni per le analisi di arricchimento
# (Reactome, GO, KEGG). Usate da three_groups_reactome.R e three_groups_go_kegg.R.
# ============================================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
})

deseq_dir <- "results/rna-seq/Exercise_GF/deseq2"
padj_cut  <- 0.05
tissues   <- c(Cb = "Cerebellum", Hc = "Hippocampus")

contrasts <- data.frame(
  label = c("Bex3KO vs WT Control", "Bex3KO EE vs WT Control", "Bex3KO EE vs Bex3KO"),
  file  = c("Bex3KO_vs_WT", "Bex3KO_EE_vs_WT", "Bex3KO_EE_vs_Bex3KO"),
  padj  = c("padj_Bex3KO_vs_WT", "padj_Bex3KO_EE_vs_WT", "padj_Bex3KO_EE_vs_Bex3KO")
)
patterns <- c("KO altered, rescued by EE", "KO altered, EE intermediate",
              "KO altered, not rescued", "EE-specific", "KO vs EE only (WT in between)")

# Nome di file a partire dal nome di un pattern
pattern_file <- function(p) gsub("[ ,()]+", "_", tolower(p))

# ID Ensembl -> ID Entrez (Reactome e KEGG lavorano con gli Entrez)
to_entrez <- function(ens) {
  m <- suppressMessages(bitr(ens, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
  unique(m$ENTREZID)
}

# Dot plot delle 20 vie/termini con p-value piu' basso (stile report del laboratorio)
enrichment_dotplot <- function(res, title, subtitle, out_png) {
  top <- head(res[order(res$pvalue), ], 20)
  top$GeneRatio_num <- sapply(strsplit(top$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  # Nomi lunghi su piu' righe, cosi' non schiacciano il grafico
  top$Description <- sapply(top$Description, function(d) paste(strwrap(d, width = 45), collapse = "\n"))
  top$Description <- factor(top$Description, levels = rev(unique(top$Description)))
  p <- ggplot(top, aes(GeneRatio_num, Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "firebrick", high = "steelblue", name = "p.adjust") +
    scale_size_continuous(range = c(3, 9)) +
    labs(title = title, subtitle = subtitle, x = "GeneRatio", y = NULL) +
    theme_bw() + theme(axis.text.y = element_text(size = 10, lineheight = 0.85),
                       plot.title.position = "plot")
  ggsave(out_png, p, width = 10, height = 9)
}
