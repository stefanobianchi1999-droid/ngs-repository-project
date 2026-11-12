# ============================================
# Exercise_GF: enrichment plot GSEA (running score) delle vie chiave, richiesti dal prof.
# Vie KEGG: Ribosome, Oxidative phosphorylation, GABAergic synapse;
# termini GO BP: cytoplasmic translation, aerobic respiration (chemical synaptic transmission non viene testato).
# Rifa' la GSEA con gli stessi parametri e la stessa classifica di three_groups_gsea.R
# (geneSet di tutte le dimensioni 15-500, pvalueCutoff = 1 per poter disegnare anche le
# vie non significative).
# Output: three_groups_<tessuto>/gsea/enrichment_plots/
# KEGG si scarica le vie dal sito KEGG: serve la connessione a internet.
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_gsea_plots.R")
# ============================================
source("pipelines/rna-seq/scripts/r/enrichment_helpers.R")
suppressPackageStartupMessages(library(enrichplot))

set.seed(1)
contrasts$stat <- c("stat_Bex3KO_vs_WT", "stat_Bex3KO_EE_vs_WT", "stat_Bex3KO_EE_vs_Bex3KO")
kegg_sets <- c(Ribosome = "mmu03010", OXPHOS = "mmu00190", GABAergic_synapse = "mmu04727")
go_sets   <- c(cytoplasmic_translation = "GO:0002181", aerobic_respiration = "GO:0009060",
               chemical_synaptic_transmission = "GO:0007268")

# Stessa classifica di three_groups_gsea.R
ranked_list <- function(tab, stat_col, to = c("ENSEMBL", "ENTREZID")) {
  to <- match.arg(to)
  d <- tab[!is.na(tab[[stat_col]]), c("gene_id", stat_col)]
  if (to == "ENTREZID") {
    m <- suppressMessages(bitr(d$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
    d <- merge(d, m, by.x = "gene_id", by.y = "ENSEMBL")
    d <- d[order(-abs(d[[stat_col]])), ]
    d <- d[!duplicated(d$ENTREZID), ]
    ids <- d$ENTREZID
  } else {
    ids <- d$gene_id
  }
  sort(setNames(d[[stat_col]], ids), decreasing = TRUE)
}

# Un enrichment plot per ogni via; il sottotitolo riporta NES e p.adjust
plot_sets <- function(gs, sets, db, title, out_dir, file_tag) {
  res <- as.data.frame(gs)
  for (nm in names(sets)) {
    r <- res[res$ID == sets[[nm]], ]
    if (nrow(r) == 0) {
      cat("  ", db, nm, title, ": non testata\n")
      next
    }
    label <- sub(" - Mus musculus \\(house mouse\\)$", "", r$Description)
    p <- gseaplot2(gs, geneSetID = r$ID,
                   title = sprintf("%s (%s) | %s\nNES = %.2f, p.adjust = %.2g",
                                   label, db, title, r$NES, r$p.adjust))
    ggsave(file.path(out_dir, paste0(db, "_", nm, "_", file_tag, ".png")), p, width = 9, height = 6)
    cat("  ", db, nm, title, ": NES", round(r$NES, 2), "p.adjust", signif(r$p.adjust, 2), "\n")
  }
}

for (t in names(tissues)) {
  cat("\n==========", tissues[[t]], "==========\n")
  in_dir  <- file.path(deseq_dir, paste0("three_groups_", t))
  out_dir <- file.path(in_dir, "gsea", "enrichment_plots")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  tab <- read.csv(file.path(in_dir, "three_groups_results.csv"))

  for (i in seq_len(nrow(contrasts))) {
    title  <- paste0(contrasts$label[i], ", ", tissues[[t]])
    rk_ens <- ranked_list(tab, contrasts$stat[i], "ENSEMBL")
    rk_ent <- ranked_list(tab, contrasts$stat[i], "ENTREZID")

    kg <- gseKEGG(rk_ent, organism = "mmu", minGSSize = 15, maxGSSize = 500,
                  pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)
    plot_sets(kg, kegg_sets, "KEGG", title, out_dir, contrasts$file[i])

    go <- gseGO(rk_ens, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP",
                minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, eps = 0,
                seed = TRUE, verbose = FALSE)
    plot_sets(go, go_sets, "GO", title, out_dir, contrasts$file[i])
  }
}
