# ============================================
# Exercise_GF: GSEA (Gene Set Enrichment Analysis) per ogni tessuto e confronto
# del modello a 3 gruppi, su GO Biological Process, KEGG e Reactome.
# A differenza dell'arricchimento sui DEG usa TUTTI i geni, ordinati per la
# statistica di Wald di DESeq2 (stat_*): positivo = piu' espresso nel primo gruppo.
# NES > 0 = via attivata nel primo gruppo del confronto, NES < 0 = via repressa.
# Output: three_groups_<tessuto>/gsea/ (tabelle + dot plot + running score di esempio)
# KEGG si scarica le vie dal sito KEGG: serve la connessione a internet.
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_gsea.R")
# ============================================
source("pipelines/rna-seq/scripts/r/enrichment_helpers.R")
suppressPackageStartupMessages({
  library(ReactomePA)
  library(enrichplot)
})

set.seed(1)   # la GSEA usa permutazioni casuali: il seme rende i risultati riproducibili
contrasts$stat <- c("stat_Bex3KO_vs_WT", "stat_Bex3KO_EE_vs_WT", "stat_Bex3KO_EE_vs_Bex3KO")

# Lista ordinata (dal piu' "su" al piu' "giu'") con nomi = ID dei geni
ranked_list <- function(tab, stat_col, to = c("ENSEMBL", "ENTREZID")) {
  to <- match.arg(to)
  d <- tab[!is.na(tab[[stat_col]]), c("gene_id", stat_col)]
  if (to == "ENTREZID") {
    m <- suppressMessages(bitr(d$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
    d <- merge(d, m, by.x = "gene_id", by.y = "ENSEMBL")
    # Se piu' geni Ensembl puntano allo stesso Entrez, tengo quello con |stat| piu' grande
    d <- d[order(-abs(d[[stat_col]])), ]
    d <- d[!duplicated(d$ENTREZID), ]
    ids <- d$ENTREZID
  } else {
    ids <- d$gene_id
  }
  v <- setNames(d[[stat_col]], ids)
  sort(v, decreasing = TRUE)
}

# Dot plot GSEA: le 10 vie piu' significative attivate (NES > 0) e le 10 represse (NES < 0)
gsea_dotplot <- function(res, title, subtitle, out_png) {
  res <- res[order(res$p.adjust), ]
  top <- rbind(head(res[res$NES > 0, ], 10), head(res[res$NES < 0, ], 10))
  top$Description <- sapply(top$Description, function(d) paste(strwrap(d, width = 45), collapse = "\n"))
  top <- top[order(top$NES), ]
  top$Description <- factor(top$Description, levels = unique(top$Description))
  top$direction <- ifelse(top$NES > 0, "Activated", "Suppressed")
  p <- ggplot(top, aes(NES, Description, size = setSize, color = p.adjust)) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_point() +
    scale_color_gradient(low = "firebrick", high = "steelblue", name = "p.adjust") +
    scale_size_continuous(range = c(3, 8), name = "Set size") +
    labs(title = title, subtitle = subtitle, x = "Normalized enrichment score (NES)", y = NULL) +
    theme_bw() + theme(axis.text.y = element_text(size = 10, lineheight = 0.85),
                       plot.title.position = "plot")
  ggsave(out_png, p, width = 10, height = 9)
}

save_gsea <- function(gs, db_label, title, out_prefix) {
  res <- if (is.null(gs)) data.frame() else as.data.frame(gs)
  if (nrow(res) > 0 && db_label == "KEGG")
    res$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", res$Description)
  write.csv(res, paste0(out_prefix, ".csv"), row.names = FALSE)
  cat("  ", db_label, title, ":", sum(res$NES > 0), "attivate,", sum(res$NES < 0),
      "represse (p.adjust < 0.05)\n")
  if (nrow(res) > 0)
    gsea_dotplot(res, title, paste0("GSEA, ", db_label,
                                    "; genes ranked by DESeq2 Wald statistic; NES > 0 = higher in ",
                                    sub(" vs .*", "", title)),
                 paste0(out_prefix, ".png"))
  invisible(gs)
}

for (t in names(tissues)) {
  cat("\n==========", tissues[[t]], "==========\n")
  in_dir  <- file.path(deseq_dir, paste0("three_groups_", t))
  out_dir <- file.path(in_dir, "gsea")
  dir.create(out_dir, showWarnings = FALSE)
  tab <- read.csv(file.path(in_dir, "three_groups_results.csv"))

  for (i in seq_len(nrow(contrasts))) {
    title  <- paste0(contrasts$label[i], " (", tissues[[t]], ")")
    prefix <- file.path(out_dir, paste0("gsea_%s_", contrasts$file[i]))
    rk_ens <- ranked_list(tab, contrasts$stat[i], "ENSEMBL")
    rk_ent <- ranked_list(tab, contrasts$stat[i], "ENTREZID")

    # GO Biological Process (termini quasi uguali uniti con simplify)
    go <- gseGO(rk_ens, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP",
                minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0,
                seed = TRUE, verbose = FALSE)
    if (!is.null(go) && nrow(as.data.frame(go)) > 0) go <- simplify(go, cutoff = 0.7)
    if (!is.null(go)) go <- setReadable(go, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
    save_gsea(go, "GO BP", title, sprintf(prefix, "go_bp"))

    # KEGG
    kg <- gseKEGG(rk_ent, organism = "mmu", minGSSize = 15, maxGSSize = 500,
                  pvalueCutoff = 0.05, eps = 0, seed = TRUE, verbose = FALSE)
    if (!is.null(kg)) kg <- setReadable(kg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    save_gsea(kg, "KEGG", title, sprintf(prefix, "kegg"))

    # Reactome (gsePathway non accetta seed/eps: fisso il seme a mano per la riproducibilita')
    set.seed(1)
    rc <- gsePathway(rk_ent, organism = "mouse", minGSSize = 15, maxGSSize = 500,
                     pvalueCutoff = 0.05, verbose = FALSE)
    if (!is.null(rc)) rc <- setReadable(rc, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    save_gsea(rc, "Reactome", title, sprintf(prefix, "reactome"))

    # Running enrichment score di esempio: la via KEGG piu' significativa del confronto
    kres <- if (is.null(kg)) data.frame() else as.data.frame(kg)
    if (nrow(kres) > 0) {
      best <- kres$ID[which.min(kres$p.adjust)]
      p <- gseaplot2(kg, geneSetID = best,
                     title = paste0(sub(" - Mus musculus \\(house mouse\\)$", "",
                                        kres$Description[kres$ID == best]), " — ", title))
      ggsave(sprintf(paste0(prefix, "_top_running_score.png"), "kegg"), p, width = 9, height = 6)
    }
  }
}
