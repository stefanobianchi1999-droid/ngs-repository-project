# ============================================
# Exercise_GF: arricchimento GO (Biological Process) e KEGG dei DEG
# del modello a 3 gruppi, per ogni tessuto, confronto e pattern.
# DEG = padj < 0.05 (su e giu' insieme); universo = tutti i geni testati nel tessuto.
# Output: three_groups_<tessuto>/go/ e three_groups_<tessuto>/kegg/ (tabelle + dot plot)
# KEGG si scarica le vie dal sito KEGG: serve la connessione a internet.
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_go_kegg.R")
# ============================================
source("pipelines/rna-seq/scripts/r/enrichment_helpers.R")

# GO Biological Process; simplify() unisce i termini GO quasi uguali (GO e' molto ridondante)
run_go <- function(genes_ens, universe_ens, title, out_prefix) {
  if (length(genes_ens) < 10) {
    cat("  GO  ", title, ": solo", length(genes_ens), "geni, saltato\n")
    return(invisible(NULL))
  }
  ego <- enrichGO(gene = genes_ens, universe = universe_ens, OrgDb = org.Mm.eg.db,
                  keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 1,
                  readable = TRUE)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) ego <- simplify(ego, cutoff = 0.7)
  res <- if (is.null(ego)) data.frame() else as.data.frame(ego)
  write.csv(res, paste0(out_prefix, ".csv"), row.names = FALSE)
  cat("  GO  ", title, ":", length(genes_ens), "geni,", nrow(res), "termini con p.adjust < 0.05\n")
  if (nrow(res) > 0)
    enrichment_dotplot(res, title, paste0("GO Biological Process, ", length(genes_ens),
                                          " DEGs (padj < 0.05); top terms by p-value"),
                       paste0(out_prefix, ".png"))
}

run_kegg <- function(genes_ens, universe_entrez, title, out_prefix) {
  genes <- to_entrez(genes_ens)
  if (length(genes) < 10) {
    cat("  KEGG", title, ": solo", length(genes), "geni, saltato\n")
    return(invisible(NULL))
  }
  ek <- enrichKEGG(gene = genes, organism = "mmu", universe = universe_entrez,
                   pvalueCutoff = 0.05, qvalueCutoff = 1)
  if (!is.null(ek)) ek <- setReadable(ek, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  res <- if (is.null(ek)) data.frame() else as.data.frame(ek)
  # I nomi KEGG finiscono con " - Mus musculus (house mouse)": lo tolgo
  if (nrow(res) > 0) res$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", res$Description)
  write.csv(res, paste0(out_prefix, ".csv"), row.names = FALSE)
  cat("  KEGG", title, ":", length(genes), "geni,", nrow(res), "vie con p.adjust < 0.05\n")
  if (nrow(res) > 0)
    enrichment_dotplot(res, title, paste0("KEGG pathways, ", length(genes),
                                          " DEGs (padj < 0.05); top pathways by p-value"),
                       paste0(out_prefix, ".png"))
}

for (t in names(tissues)) {
  cat("\n==========", tissues[[t]], "==========\n")
  in_dir <- file.path(deseq_dir, paste0("three_groups_", t))
  go_dir <- file.path(in_dir, "go");   dir.create(go_dir, showWarnings = FALSE)
  kg_dir <- file.path(in_dir, "kegg"); dir.create(kg_dir, showWarnings = FALSE)
  tab <- read.csv(file.path(in_dir, "three_groups_results.csv"))
  universe_ens    <- tab$gene_id[!is.na(tab$LRT_padj)]
  universe_entrez <- to_entrez(universe_ens)

  # Un arricchimento per confronto e uno per pattern
  sets <- list()
  for (i in seq_len(nrow(contrasts))) {
    col <- contrasts$padj[i]
    sets[[contrasts$file[i]]] <- list(
      genes = tab$gene_id[!is.na(tab[[col]]) & tab[[col]] < padj_cut],
      title = paste0(contrasts$label[i], " (", tissues[[t]], ")"))
  }
  for (p in patterns) {
    sets[[paste0("pattern_", pattern_file(p))]] <- list(
      genes = tab$gene_id[tab$pattern == p],
      title = paste0(p, " (", tissues[[t]], ")"))
  }

  for (nm in names(sets)) {
    s <- sets[[nm]]
    run_go(s$genes, universe_ens, s$title, file.path(go_dir, paste0("go_bp_", nm)))
    run_kegg(s$genes, universe_entrez, s$title, file.path(kg_dir, paste0("kegg_", nm)))
  }
}
