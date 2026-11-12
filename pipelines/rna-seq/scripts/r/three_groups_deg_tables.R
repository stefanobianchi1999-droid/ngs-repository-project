# ============================================
# Exercise_GF: tabelle complete dei DEG in Excel (una cartella di lavoro)
#   - README con la spiegazione delle colonne
#   - un foglio per tessuto x confronto: TUTTI i geni testati con le statistiche
#     DESeq2, ordinati per padj, con una colonna che segnala i DEG (padj < 0.05)
#   - un foglio per tessuto con tutti i confronti, il test LRT e il pattern
# Legge le tabelle create da three_groups_exercise_gf.R (non rilancia DESeq2).
# Output: results/rna-seq/Exercise_GF/deliverables/Bex3KO_DEG_tables.xlsx
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_deg_tables.R")
# ============================================
library(openxlsx)

deseq_dir <- "results/rna-seq/Exercise_GF/deseq2"
out_dir   <- "results/rna-seq/Exercise_GF/deliverables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
padj_cut  <- 0.05
tissues   <- c(Cb = "Cerebellum", Hc = "Hippocampus")

contrasts <- data.frame(
  label = c("Bex3KO vs WT Control", "Bex3KO EE vs WT Control", "Bex3KO EE vs Bex3KO"),
  sheet = c("KO_vs_WT", "KO-EE_vs_WT", "KO-EE_vs_KO"),
  suffix = c("Bex3KO_vs_WT", "Bex3KO_EE_vs_WT", "Bex3KO_EE_vs_Bex3KO")
)

wb <- createWorkbook()
header_style <- createStyle(textDecoration = "bold", fgFill = "#16233B", fontColour = "#FFFFFF",
                            halign = "center", wrapText = TRUE)
sig_style    <- createStyle(fgFill = "#FCE8E3")

# ---------- README ----------
readme <- data.frame(Item = c(
  "Study", "Model", "Groups", "Tissues", "Analysis",
  "DEG definition", "log2FoldChange", "",
  "Column: gene_id", "Column: gene_name", "Column: chromosome", "Column: baseMean",
  "Column: log2FoldChange", "Column: lfcSE", "Column: stat", "Column: pvalue", "Column: padj",
  "Column: DEG", "Column: direction",
  "Sheets <tissue>_<comparison>", "Sheets <tissue>_all_comparisons", "Column: LRT_padj", "Column: pattern"),
  Description = c(
  "Bex3 knockout mice and environmental enrichment (EE), bulk RNA-seq",
  "Mus musculus, GRCm38 (mm10) transcriptome, Ensembl release 102",
  "WT Control, Bex3KO, Bex3KO EE (n = 3 per group and tissue)",
  "Cerebellum (Cb), Hippocampus (Hc), analysed separately",
  "Salmon quantification, tximport gene-level summarisation, DESeq2 model ~ group per tissue (pairwise Wald tests)",
  "Adjusted p-value (Benjamini-Hochberg) < 0.05",
  "Positive = higher in the first group of the comparison (e.g. Bex3KO in 'Bex3KO vs WT Control')",
  "",
  "Ensembl gene ID", "Gene symbol", "Chromosome (Y-chromosome genes were excluded from the analysis)",
  "Mean of normalised counts across the samples of the tissue",
  "log2 fold change (first group / second group)", "Standard error of the log2 fold change",
  "Wald statistic (log2FoldChange / lfcSE); used to rank genes for GSEA",
  "Raw p-value", "Benjamini-Hochberg adjusted p-value (NA = filtered by DESeq2 independent filtering or outlier)",
  "TRUE if padj < 0.05 (highlighted rows)", "Up / Down in the first group, for DEGs only",
  "All tested genes for one comparison, sorted by padj",
  "All comparisons side by side, plus the LRT test and the expression pattern",
  "Likelihood ratio test: gene changes among at least two of the three groups",
  "Expression pattern across the three groups (e.g. 'KO altered, rescued by EE'), for LRT-significant genes")
)
addWorksheet(wb, "README")
writeData(wb, "README", readme, headerStyle = header_style)
setColWidths(wb, "README", cols = 1:2, widths = c(32, 120))

# ---------- Un foglio per tessuto x confronto ----------
summary_rows <- list()
for (t in names(tissues)) {
  tab <- read.csv(file.path(deseq_dir, paste0("three_groups_", t), "three_groups_results.csv"))

  for (i in seq_len(nrow(contrasts))) {
    s <- contrasts$suffix[i]
    d <- data.frame(
      gene_id        = tab$gene_id,
      gene_name      = tab$gene_name,
      chromosome     = tab$chrom,
      baseMean       = tab$baseMean,
      log2FoldChange = tab[[paste0("log2FC_", s)]],
      lfcSE          = tab[[paste0("lfcSE_", s)]],
      stat           = tab[[paste0("stat_", s)]],
      pvalue         = tab[[paste0("pvalue_", s)]],
      padj           = tab[[paste0("padj_", s)]]
    )
    d$DEG <- !is.na(d$padj) & d$padj < padj_cut
    d$direction <- ifelse(d$DEG, ifelse(d$log2FoldChange > 0, "Up", "Down"), "")
    d <- d[order(d$padj, na.last = TRUE), ]

    sheet <- paste0(t, "_", contrasts$sheet[i])
    addWorksheet(wb, sheet)
    writeData(wb, sheet, paste0(tissues[[t]], ": ", contrasts$label[i],
                                " (log2FC > 0 = higher in ", sub(" vs .*", "", contrasts$label[i]),
                                "); DEGs: padj < 0.05, highlighted"), startRow = 1)
    writeData(wb, sheet, d, startRow = 2, headerStyle = header_style)
    n_deg <- sum(d$DEG)
    if (n_deg > 0) addStyle(wb, sheet, sig_style, rows = 2 + seq_len(n_deg), cols = 1:ncol(d),
                            gridExpand = TRUE)
    addStyle(wb, sheet, createStyle(numFmt = "0.00"), rows = 3:(nrow(d) + 2), cols = 4:7,
             gridExpand = TRUE, stack = TRUE)
    addStyle(wb, sheet, createStyle(numFmt = "0.00E+00"), rows = 3:(nrow(d) + 2), cols = 8:9,
             gridExpand = TRUE, stack = TRUE)
    freezePane(wb, sheet, firstActiveRow = 3, firstActiveCol = 3)
    addFilter(wb, sheet, rows = 2, cols = 1:ncol(d))
    setColWidths(wb, sheet, cols = 1:ncol(d), widths = c(22, 16, 11, 11, 15, 9, 9, 11, 11, 8, 10))

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Tissue = tissues[[t]], Comparison = contrasts$label[i], Sheet = sheet,
      Genes_tested = sum(!is.na(d$padj)), DEGs = n_deg,
      Up = sum(d$direction == "Up"), Down = sum(d$direction == "Down"))
  }

  # ---------- Foglio con tutti i confronti + LRT + pattern ----------
  all <- tab
  names(all)[names(all) == "chrom"] <- "chromosome"
  sheet <- paste0(t, "_all_comparisons")
  addWorksheet(wb, sheet)
  writeData(wb, sheet, all, headerStyle = header_style)
  freezePane(wb, sheet, firstActiveRow = 2, firstActiveCol = 3)
  addFilter(wb, sheet, rows = 1, cols = 1:ncol(all))
  setColWidths(wb, sheet, cols = 1:ncol(all), widths = "auto")
}

# ---------- Riepilogo dei numeri, subito dopo il README ----------
summ <- do.call(rbind, summary_rows)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summ, headerStyle = header_style)
setColWidths(wb, "Summary", cols = 1:ncol(summ), widths = c(14, 26, 18, 14, 8, 8, 8))
worksheetOrder(wb) <- c(1, length(names(wb)), 2:(length(names(wb)) - 1))

out_file <- file.path(out_dir, "Bex3KO_DEG_tables.xlsx")
saveWorkbook(wb, out_file, overwrite = TRUE)
print(summ, row.names = FALSE)
cat("Saved:", out_file, "\n")
