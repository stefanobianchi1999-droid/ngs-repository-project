# ============================================
# Exercise_GF: tabelle supplementari S3-S5 in Excel (stesso stile di Bex3KO_DEG_tables.xlsx)
#   - S3 Bex3KO_ORA_tables.xlsx: GO BP, KEGG, Reactome dei DEG di ogni confronto e
#     dei geni di ogni pattern (un foglio per tessuto x database)
#   - S4 Bex3KO_GSEA_tables.xlsx: GSEA GO BP, KEGG, Reactome di ogni confronto
#   - S5 Bex3KO_normalized_counts.xlsx: conte normalizzate DESeq2, un foglio per tessuto
# Legge i CSV gia' creati da three_groups_go_kegg.R, three_groups_reactome.R,
# three_groups_gsea.R e three_groups_exercise_gf.R (non rilancia le analisi).
# Output: results/rna-seq/Exercise_GF/deliverables/
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_supp_tables.R")
# ============================================
library(openxlsx)

deseq_dir <- "results/rna-seq/Exercise_GF/deseq2"
out_dir   <- "results/rna-seq/Exercise_GF/deliverables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
padj_cut  <- 0.05
tissues   <- c(Cb = "Cerebellum", Hc = "Hippocampus")

contrasts <- data.frame(
  label = c("Bex3KO vs WT Control", "Bex3KO EE vs WT Control", "Bex3KO EE vs Bex3KO"),
  file  = c("Bex3KO_vs_WT", "Bex3KO_EE_vs_WT", "Bex3KO_EE_vs_Bex3KO")
)
patterns <- c("KO altered, rescued by EE", "KO altered, EE intermediate",
              "KO altered, not rescued", "EE-specific", "KO vs EE only (WT in between)")
pattern_file <- function(p) gsub("[ ,()]+", "_", tolower(p))

# Database: cartella, prefisso dei file ORA, prefisso dei file GSEA
dbs <- data.frame(
  name   = c("GO_BP", "KEGG", "Reactome"),
  folder = c("go", "kegg", "reactome"),
  ora    = c("go_bp_", "kegg_", "reactome_"),
  gsea   = c("gsea_go_bp_", "gsea_kegg_", "gsea_reactome_")
)

header_style <- createStyle(textDecoration = "bold", fgFill = "#16233B", fontColour = "#FFFFFF",
                            halign = "center", wrapText = TRUE)
num2 <- createStyle(numFmt = "0.00")
sci  <- createStyle(numFmt = "0.00E+00")

# CSV -> data.frame; NULL se il file non esiste (insieme saltato perche' < 10 geni)
read_res <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, check.names = FALSE)
}

readme_sheet <- function(wb, items) {
  addWorksheet(wb, "README")
  writeData(wb, "README", data.frame(Item = names(items), Description = unname(items)),
            headerStyle = header_style)
  setColWidths(wb, "README", cols = 1:2, widths = c(32, 120))
}

# Scrive un foglio di risultati con intestazione, filtri e formati numerici
result_sheet <- function(wb, sheet, title, d, cols_num2, cols_sci, widths) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, title, startRow = 1)
  writeData(wb, sheet, d, startRow = 2, headerStyle = header_style)
  if (nrow(d) > 0) {
    rows <- 3:(nrow(d) + 2)
    addStyle(wb, sheet, num2, rows = rows, cols = which(names(d) %in% cols_num2), gridExpand = TRUE)
    addStyle(wb, sheet, sci,  rows = rows, cols = which(names(d) %in% cols_sci),  gridExpand = TRUE)
  }
  freezePane(wb, sheet, firstActiveRow = 3, firstActiveCol = 4)
  addFilter(wb, sheet, rows = 2, cols = 1:ncol(d))
  setColWidths(wb, sheet, cols = 1:ncol(d), widths = widths)
}

# Numero di geni di input di ogni insieme (DEG di un confronto o geni di un pattern)
input_sets <- function(tab) {
  s <- list()
  for (i in seq_len(nrow(contrasts))) {
    col <- paste0("padj_", contrasts$file[i])
    s[[contrasts$label[i]]] <- list(file = contrasts$file[i], type = "DEGs of comparison",
                                    n = sum(!is.na(tab[[col]]) & tab[[col]] < padj_cut))
  }
  for (p in patterns)
    s[[p]] <- list(file = paste0("pattern_", pattern_file(p)), type = "Expression pattern",
                   n = sum(tab$pattern == p, na.rm = TRUE))
  s
}

# ==================== S3: over-representation analysis ====================
wb <- createWorkbook()
readme_sheet(wb, c(
  "Table" = "Supplementary Table S3. Over-representation analysis (ORA)",
  "Input gene sets" = "DEGs (padj < 0.05, up and down together) of each comparison, and genes of each expression pattern",
  "Background" = "All genes tested in the tissue",
  "GO_BP" = "clusterProfiler enrichGO, Biological Process; redundant terms removed with simplify (cutoff 0.7)",
  "KEGG" = "clusterProfiler enrichKEGG, organism mmu (KEGG database accessed online)",
  "Reactome" = "ReactomePA enrichPathway, organism mouse",
  "Significance" = "Benjamini-Hochberg adjusted p-value < 0.05; only significant terms are listed",
  " " = "",
  "Sheet Summary" = "Number of input genes and significant terms for every tissue, gene set and database",
  "Sheets <tissue>_<database>" = "All significant terms of one tissue and database; filter column Gene_set to see one comparison or pattern",
  "Column: Gene_set" = "Comparison or expression pattern used as input",
  "Column: GeneRatio" = "Input genes in the term / input genes annotated in the database",
  "Column: BgRatio" = "Background genes in the term / background genes annotated in the database",
  "Column: FoldEnrichment" = "GeneRatio / BgRatio",
  "Column: pvalue, p.adjust, qvalue" = "Hypergeometric test p-value, Benjamini-Hochberg adjusted p-value, q-value",
  "Column: Count" = "Number of input genes in the term",
  "Column: geneID" = "Input genes in the term (gene symbols)",
  "Note" = "KEGG disease maps (e.g. Parkinson disease, Coronavirus disease) share many genes with Oxidative phosphorylation and Ribosome: check column geneID before interpreting them"
))

summary_rows <- list()
for (t in names(tissues)) {
  tab  <- read.csv(file.path(deseq_dir, paste0("three_groups_", t), "three_groups_results.csv"))
  sets <- input_sets(tab)
  for (k in seq_len(nrow(dbs))) {
    parts <- list()
    for (nm in names(sets)) {
      s   <- sets[[nm]]
      res <- read_res(file.path(deseq_dir, paste0("three_groups_", t), dbs$folder[k],
                                paste0(dbs$ora[k], s$file, ".csv")))
      status <- if (is.null(res)) "not run (< 10 genes)" else as.character(nrow(res))
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Tissue = tissues[[t]], Database = dbs$name[k], Gene_set = nm, Type = s$type,
        Input_genes = s$n, Significant_terms = status)
      if (!is.null(res) && nrow(res) > 0) parts[[nm]] <- cbind(Gene_set = nm, res)
    }
    keep <- c("Gene_set", "category", "subcategory", "ID", "Description", "GeneRatio", "BgRatio",
              "FoldEnrichment", "pvalue", "p.adjust", "qvalue", "Count", "geneID")
    d <- do.call(rbind, lapply(parts, function(x) x[, intersect(keep, names(x))]))
    rownames(d) <- NULL
    widths <- ifelse(names(d) %in% c("Description", "geneID"), 50,
              ifelse(names(d) %in% c("Gene_set", "category", "subcategory"), 28, 12))
    result_sheet(wb, paste0(t, "_", dbs$name[k]),
                 paste0(tissues[[t]], ": ", dbs$name[k], " over-representation (p.adjust < 0.05)"),
                 d, "FoldEnrichment", c("pvalue", "p.adjust", "qvalue"), widths)
  }
}
summ <- do.call(rbind, summary_rows)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summ, headerStyle = header_style)
setColWidths(wb, "Summary", cols = 1:ncol(summ), widths = c(14, 11, 32, 20, 12, 20))
worksheetOrder(wb) <- c(1, length(names(wb)), 2:(length(names(wb)) - 1))
saveWorkbook(wb, file.path(out_dir, "Bex3KO_ORA_tables.xlsx"), overwrite = TRUE)
cat("Saved:", file.path(out_dir, "Bex3KO_ORA_tables.xlsx"), "\n")

# ==================== S4: GSEA ====================
wb <- createWorkbook()
readme_sheet(wb, c(
  "Table" = "Supplementary Table S4. Gene set enrichment analysis (GSEA)",
  "Ranking" = "All tested genes ranked by the DESeq2 Wald statistic of each comparison",
  "Methods" = "clusterProfiler gseGO (Biological Process), gseKEGG (organism mmu), ReactomePA gsePathway; gene-set size 15-500",
  "Significance" = "Benjamini-Hochberg adjusted p-value < 0.05; only significant gene sets are listed",
  " " = "",
  "Sheet Summary" = "Number of significant gene sets, higher and lower in the first group, for every tissue, comparison and database",
  "Sheets <tissue>_<database>" = "All significant gene sets of one tissue and database; filter column Comparison to see one comparison",
  "Column: Direction" = "Higher / Lower in the first group of the comparison (sign of NES)",
  "Column: setSize" = "Number of genes of the gene set among the ranked genes",
  "Column: NES" = "Normalised enrichment score; positive = higher in the first group (e.g. Bex3KO in 'Bex3KO vs WT Control')",
  "Column: pvalue, p.adjust, qvalue" = "Permutation p-value, Benjamini-Hochberg adjusted p-value, q-value",
  "Column: leading_edge" = "tags: % of the gene set in the leading edge; list: % of the ranked list before the peak; signal: combined",
  "Column: core_enrichment" = "Leading-edge genes (gene symbols): the genes that drive the enrichment"
))

summary_rows <- list()
for (t in names(tissues)) {
  for (k in seq_len(nrow(dbs))) {
    parts <- list()
    for (i in seq_len(nrow(contrasts))) {
      res <- read_res(file.path(deseq_dir, paste0("three_groups_", t), "gsea",
                                paste0(dbs$gsea[k], contrasts$file[i], ".csv")))
      if (is.null(res)) res <- data.frame(NES = numeric(0))
      first <- sub(" vs .*", "", contrasts$label[i])
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Tissue = tissues[[t]], Database = dbs$name[k], Comparison = contrasts$label[i],
        Significant_sets = nrow(res),
        Higher_in_first = sum(res$NES > 0), Lower_in_first = sum(res$NES < 0))
      if (nrow(res) > 0) {
        res$Direction <- ifelse(res$NES > 0, paste("Higher in", first), paste("Lower in", first))
        parts[[contrasts$label[i]]] <- cbind(Comparison = contrasts$label[i], res)
      }
    }
    keep <- c("Comparison", "ID", "Description", "Direction", "setSize", "NES",
              "pvalue", "p.adjust", "qvalue", "leading_edge", "core_enrichment")
    d <- do.call(rbind, lapply(parts, function(x) x[, keep]))
    rownames(d) <- NULL
    widths <- c(24, 16, 50, 22, 9, 9, 11, 11, 11, 28, 60)
    result_sheet(wb, paste0(t, "_", dbs$name[k]),
                 paste0(tissues[[t]], ": GSEA of ", dbs$name[k], " gene sets (p.adjust < 0.05)"),
                 d, "NES", c("pvalue", "p.adjust", "qvalue"), widths)
  }
}
summ <- do.call(rbind, summary_rows)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summ, headerStyle = header_style)
setColWidths(wb, "Summary", cols = 1:ncol(summ), widths = c(14, 11, 26, 16, 16, 16))
worksheetOrder(wb) <- c(1, length(names(wb)), 2:(length(names(wb)) - 1))
saveWorkbook(wb, file.path(out_dir, "Bex3KO_GSEA_tables.xlsx"), overwrite = TRUE)
cat("Saved:", file.path(out_dir, "Bex3KO_GSEA_tables.xlsx"), "\n")

# ==================== S5: conte normalizzate ====================
# Tabella dei campioni come in exercise_gf_setup.R; sesso dedotto dai geni del cromosoma Y
samples <- data.frame(
  np     = paste0("NP", 1:18),
  group  = rep(c("Bex3KO", "Bex3KO", "Bex3KO_EE", "Bex3KO_EE", "WT_Control", "WT_Control"), each = 3),
  tissue = rep(c("Hc", "Cb", "Hc", "Cb", "Hc", "Cb"), each = 3),
  rep    = rep(1:3, times = 6)
)
samples$sample <- paste(samples$group, samples$tissue, samples$rep, sep = "_")
samples$inferred_sex <- ifelse(samples$group == "WT_Control", "F",
                        ifelse(samples$group == "Bex3KO_EE", "M",
                        ifelse(samples$rep == 2, "F", "M")))

wb <- createWorkbook()
readme_sheet(wb, c(
  "Table" = "Supplementary Table S5. DESeq2-normalised counts",
  "Values" = "Counts divided by the DESeq2 size factor of each sample (median-of-ratios), estimated separately for each tissue",
  "Genes" = "All genes analysed in the tissue (at least 10 counts in at least 3 samples; Y-chromosome genes excluded)",
  "Sheet Samples" = "Sample names used in the column headers, sequencing ID (NP), group, tissue, replicate and sex inferred from Y-chromosome genes",
  "Sheets Cb, Hc" = "One row per gene: Ensembl ID, gene symbol, normalised counts per sample",
  "Raw data" = "Raw FASTQ files and raw counts: GEO [accession GSE...]"
))
addWorksheet(wb, "Samples")
writeData(wb, "Samples", samples[, c("sample", "np", "group", "tissue", "rep", "inferred_sex")],
          headerStyle = header_style)
setColWidths(wb, "Samples", cols = 1:6, widths = c(20, 8, 13, 8, 6, 13))

for (t in names(tissues)) {
  nc <- read.csv(file.path(deseq_dir, paste0("three_groups_", t), "normalized_counts.csv"),
                 check.names = FALSE)
  names(nc)[1] <- "gene_id"
  d <- nc[, c("gene_id", "gene_name", setdiff(names(nc), c("gene_id", "gene_name")))]
  addWorksheet(wb, t)
  writeData(wb, t, paste0(tissues[[t]], ": DESeq2-normalised counts"), startRow = 1)
  writeData(wb, t, d, startRow = 2, headerStyle = header_style)
  addStyle(wb, t, createStyle(numFmt = "0.0"), rows = 3:(nrow(d) + 2), cols = 3:ncol(d),
           gridExpand = TRUE)
  freezePane(wb, t, firstActiveRow = 3, firstActiveCol = 3)
  addFilter(wb, t, rows = 2, cols = 1:ncol(d))
  setColWidths(wb, t, cols = 1:ncol(d), widths = c(22, 16, rep(17, ncol(d) - 2)))
}
saveWorkbook(wb, file.path(out_dir, "Bex3KO_normalized_counts.xlsx"), overwrite = TRUE)
cat("Saved:", file.path(out_dir, "Bex3KO_normalized_counts.xlsx"), "\n")
