# ============================================
# Exercise_GF: arricchimento Reactome dei DEG del modello a 3 gruppi
# Per ogni tessuto e confronto: DEG (padj < 0.05) -> ID Entrez -> enrichPathway
# (universo = tutti i geni testati) -> tabella + dot plot (stile report del laboratorio).
# Fa lo stesso anche per i geni di ogni pattern (es. "KO altered, rescued by EE").
# Legge le tabelle create da three_groups_exercise_gf.R (non rilancia DESeq2).
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_reactome.R")
# ============================================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
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

# ID Ensembl -> ID Entrez (Reactome lavora con gli Entrez)
to_entrez <- function(ens) {
  m <- suppressMessages(bitr(ens, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
  unique(m$ENTREZID)
}

# Arricchimento + tabella + dot plot; restituisce il numero di vie significative
run_reactome <- function(genes_ens, universe_entrez, title, out_prefix) {
  genes <- to_entrez(genes_ens)
  if (length(genes) < 10) {
    cat("  ", title, ": solo", length(genes), "geni, arricchimento saltato\n")
    return(0)
  }
  er <- enrichPathway(gene = genes, organism = "mouse", universe = universe_entrez,
                      pvalueCutoff = 0.05, qvalueCutoff = 1, readable = TRUE)
  res <- if (is.null(er)) data.frame() else as.data.frame(er)
  write.csv(res, paste0(out_prefix, ".csv"), row.names = FALSE)
  # enrichPathway applica pvalueCutoff al p-value corretto (BH): qui restano le vie con p.adjust < 0.05
  n_sig <- sum(res$p.adjust < 0.05)
  cat("  ", title, ":", length(genes), "geni,", n_sig, "vie con p.adjust < 0.05\n")
  if (nrow(res) == 0) return(0)

  # Dot plot delle 20 vie con p-value piu' basso (come nei report del laboratorio)
  top <- head(res[order(res$pvalue), ], 20)
  top$GeneRatio_num <- sapply(strsplit(top$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  # Nomi lunghi delle vie su piu' righe, cosi' non schiacciano il grafico
  top$Description <- sapply(top$Description, function(d) paste(strwrap(d, width = 45), collapse = "
"))
  top$Description <- factor(top$Description, levels = rev(unique(top$Description)))
  p <- ggplot(top, aes(GeneRatio_num, Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "firebrick", high = "steelblue", name = "p.adjust") +
    scale_size_continuous(range = c(3, 9)) +
    labs(title = title, subtitle = paste0("Reactome enrichment, ", length(genes),
                                          " DEGs (padj < 0.05); top pathways by p-value"),
         x = "GeneRatio", y = NULL) +
    theme_bw() + theme(axis.text.y = element_text(size = 10, lineheight = 0.85),
                       plot.title.position = "plot")
  ggsave(paste0(out_prefix, ".png"), p, width = 10, height = 9)
  n_sig
}

for (t in names(tissues)) {
  cat("\n==========", tissues[[t]], "==========\n")
  in_dir  <- file.path(deseq_dir, paste0("three_groups_", t))
  out_dir <- file.path(in_dir, "reactome")
  dir.create(out_dir, showWarnings = FALSE)
  tab <- read.csv(file.path(in_dir, "three_groups_results.csv"))
  universe <- to_entrez(tab$gene_id[!is.na(tab$LRT_padj)])

  # Un arricchimento per confronto (tutti i DEG, su e giu' insieme)
  for (i in seq_len(nrow(contrasts))) {
    col <- contrasts$padj[i]
    deg <- tab$gene_id[!is.na(tab[[col]]) & tab[[col]] < padj_cut]
    run_reactome(deg, universe, paste0(contrasts$label[i], " (", tissues[[t]], ")"),
                 file.path(out_dir, paste0("reactome_", contrasts$file[i])))
  }
  # Un arricchimento per pattern
  for (p in patterns) {
    deg <- tab$gene_id[tab$pattern == p]
    run_reactome(deg, universe, paste0(p, " (", tissues[[t]], ")"),
                 file.path(out_dir, paste0("reactome_pattern_", gsub("[ ,()]+", "_", tolower(p)))))
  }
}
