# ============================================
# Exercise_GF: pacchetto per la sottomissione a GEO, da consegnare al laboratorio
#   - processed/: conte grezze e TPM (tutti i geni, tximport), conte normalizzate DESeq2,
#     quant.sf di Salmon rinominati con il nome del campione (compressi)
#   - GEO_metadata_draft.xlsx: metadati nella struttura del modello GEO per sequencing
#     (STUDY, SAMPLES, PROTOCOLS, DATA PROCESSING, MD5); in giallo i campi che deve
#     completare il laboratorio (FASTQ grezzi, kit, strumento, animali, protocollo EE)
# I FASTQ grezzi non sono inclusi: li carica il laboratorio.
# Output: results/rna-seq/Exercise_GF/deliverables/GEO_package/
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/geo_package.R")
# ============================================
source("pipelines/rna-seq/scripts/r/exercise_gf_setup.R")
library(openxlsx)

out_dir  <- "results/rna-seq/Exercise_GF/deliverables/GEO_package"
proc_dir <- file.path(out_dir, "processed")
dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)
tissues  <- c(Cb = "Cerebellum", Hc = "Hippocampus")

# Sesso dedotto dai geni del cromosoma Y (come in three_groups_supp_tables.R)
samples$sex <- ifelse(samples$group == "WT_Control", "female",
               ifelse(samples$group == "Bex3KO_EE", "male",
               ifelse(samples$rep == 2, "female", "male")))

# ---------- Dati processati ----------
proc_files <- character(0)
for (t in names(tissues)) {
  keep <- samples$sample[samples$tissue == t]
  for (what in c("raw_counts", "TPM")) {
    m <- if (what == "raw_counts") round(txi_all$counts[, keep]) else txi_all$abundance[, keep]
    d <- data.frame(gene_id = rownames(m), gene_name = symbol_of(rownames(m)), m, check.names = FALSE)
    f <- file.path(proc_dir, paste0("Bex3KO_", what, "_", t, ".csv.gz"))
    write.csv(d, gzfile(f), row.names = FALSE)
    proc_files <- c(proc_files, f)
  }
  f <- file.path(proc_dir, paste0("Bex3KO_normalized_counts_", t, ".csv.gz"))
  write.csv(read.csv(file.path(deseq_dir, paste0("three_groups_", t), "normalized_counts.csv"),
                     check.names = FALSE), gzfile(f), row.names = FALSE)
  proc_files <- c(proc_files, f)
}
quant_files <- file.path(proc_dir, paste0(samples$sample, "_quant.sf.gz"))
for (i in seq_len(nrow(samples))) {
  con <- gzfile(quant_files[i], "w")
  writeLines(readLines(file.path(salmon_dir, samples$np[i], "quant.sf")), con)
  close(con)
}

# ---------- Foglio dei metadati ----------
todo <- function(x) paste0("[TO COMPLETE: ", x, "]")
wb <- createWorkbook()
header_style <- createStyle(textDecoration = "bold", fgFill = "#16233B", fontColour = "#FFFFFF", wrapText = TRUE)
todo_style   <- createStyle(fgFill = "#FFF2A8")
wrap_style   <- createStyle(wrapText = TRUE, valign = "top")

# Scrive un foglio e colora in giallo le celle ancora da completare
add_sheet <- function(name, d, widths) {
  addWorksheet(wb, name)
  writeData(wb, name, d, headerStyle = header_style)
  addStyle(wb, name, wrap_style, rows = 2:(nrow(d) + 1), cols = 1:ncol(d), gridExpand = TRUE)
  for (j in seq_len(ncol(d))) {
    r <- which(grepl("^\\[TO COMPLETE", d[[j]]))
    if (length(r)) addStyle(wb, name, todo_style, rows = r + 1, cols = j, stack = TRUE)
  }
  setColWidths(wb, name, cols = seq_len(ncol(d)), widths = widths)
  freezePane(wb, name, firstActiveRow = 2)
}

add_sheet("README", data.frame(Item = c(
  "Purpose", "What is included", "What the lab must add", "Yellow cells", "Template",
  "Raw files", "Sex", "Contact"),
  Description = c(
  "Draft of the GEO submission for the Bex3KO / environmental enrichment RNA-seq study (18 samples).",
  "Processed data in folder 'processed' (raw counts, TPM and DESeq2-normalised counts per tissue; Salmon quant.sf per sample) with MD5 checksums in sheet MD5.",
  "Raw FASTQ files (untrimmed, as delivered by HaploX) for all 18 samples, their file names and MD5 checksums, and the fields marked in yellow.",
  "Information that is not available from the data (animals, EE protocol, RNA extraction, library kit, sequencer, contributors).",
  "GEO updates its metadata template: download the current 'high-throughput sequencing' template from the GEO submission page and copy these fields into it.",
  "Paired-end reads: two FASTQ files per sample (R1, R2). GEO forwards them to SRA.",
  "Inferred from Y-chromosome gene expression (Ddx3y, Eif2s3y, Kdm5d, Uty); please confirm with the animal records.",
  todo("submitter name and e-mail; the GEO record is normally submitted by the lab / principal investigator"))),
  c(24, 120))

add_sheet("STUDY", data.frame(Field = c(
  "title", "summary", "overall design", "contributor", "contributor", "supplementary file", "supplementary file",
  "supplementary file", "supplementary file", "supplementary file", "supplementary file"),
  Value = c(
  "Transcriptomic effects of Bex3 deletion and environmental enrichment in mouse hippocampus and cerebellum",
  todo("2-4 sentences: background on Bex3, aim of the study, main finding. Draft: Bulk RNA-seq of hippocampus and cerebellum from wild-type, Bex3 knockout and Bex3 knockout mice housed in an enriched environment, to identify genes altered by Bex3 loss and test whether environmental enrichment reverses these changes."),
  "Hippocampus and cerebellum were collected from WT Control, Bex3KO and Bex3KO mice housed in an enriched environment (Bex3KO EE), n = 3 animals per group (18 samples). Paired-end, stranded RNA-seq libraries were sequenced and analysed separately for each tissue.",
  todo("Surname,Initial (one row per author)"), todo("Surname,Initial"),
  basename(proc_files))),
  c(22, 120))

geo_samples <- data.frame(
  `library name`        = samples$sample,
  title                 = paste0(tissues[samples$tissue], ", ", gsub("_", " ", samples$group), ", replicate ", samples$rep),
  organism              = "Mus musculus",
  tissue                = unname(tissues[samples$tissue]),
  `cell line`           = "",
  `cell type`           = "",
  genotype              = ifelse(samples$group == "WT_Control", "wild type", "Bex3 knockout"),
  treatment             = ifelse(samples$group == "Bex3KO_EE", "environmental enrichment", "standard housing"),
  sex                   = samples$sex,
  age                   = todo("age at sacrifice"),
  strain                = todo("background strain"),
  molecule              = todo("polyA RNA or total RNA"),
  `single or paired-end` = "paired-end",
  `instrument model`    = todo("e.g. Illumina NovaSeq 6000"),
  description           = paste0("sequencing ID ", samples$np),
  `processed data file` = basename(quant_files),
  `raw file 1`          = todo(paste0(samples$np, " R1 FASTQ file name")),
  `raw file 2`          = todo(paste0(samples$np, " R2 FASTQ file name")),
  check.names = FALSE)
add_sheet("SAMPLES", geo_samples, c(18, 36, 13, 12, 9, 9, 14, 22, 8, 20, 20, 10, 12, 26, 18, 30, 30, 30))

add_sheet("PROTOCOLS", data.frame(Field = c(
  "growth protocol", "treatment protocol", "extract protocol", "library construction protocol",
  "library strategy"),
  Value = c(
  todo("animals: strain, generation of Bex3 KO, housing, age, sex"),
  todo("environmental enrichment: cage type, objects, duration, start age"),
  todo("dissection and RNA extraction method, RNA quality (RIN)"),
  todo("library kit, poly-A selection, stranded protocol (libraries detected as stranded, reverse; ISR), read length 2 x 150 bp, provider HaploX"),
  "RNA-Seq")),
  c(30, 120))

add_sheet("DATA PROCESSING", data.frame(Step = c(
  "data processing step", "data processing step", "data processing step", "data processing step",
  "data processing step", "data processing step", "genome build/assembly", "processed data files format and content"),
  Value = c(
  "Read quality was assessed with FastQC v0.13.0; adapters and low-quality bases were removed with fastp v1.4.0 (--detect_adapter_for_pe).",
  "Trimmed reads were quantified with Salmon v2.8.0 in mapping-based mode against the GRCm38 cDNA transcriptome (Ensembl release 102), with automatic library-type detection (-l A; ISR) and --gcBias.",
  "Transcript estimates were summarised to gene level with tximport v1.40.0 (R v4.6.1) using a transcript-to-gene table from the Ensembl 102 GTF.",
  "Differential expression was tested with DESeq2 v1.52.0 separately for each tissue (design ~ group; likelihood ratio test and pairwise Wald tests); genes with at least 10 counts in at least 3 samples were kept and Y-chromosome genes removed.",
  "Normalised counts are DESeq2 median-of-ratios normalised counts of the genes analysed in each tissue.",
  "Sex of each animal was inferred from Y-chromosome gene expression (Ddx3y, Eif2s3y, Kdm5d, Uty).",
  "GRCm38 (mm10), Ensembl release 102",
  "*_quant.sf.gz: Salmon transcript-level output per sample (TPM and estimated reads). Bex3KO_raw_counts_<tissue>.csv.gz: gene-level estimated counts from tximport, rounded, all genes. Bex3KO_TPM_<tissue>.csv.gz: gene-level TPM. Bex3KO_normalized_counts_<tissue>.csv.gz: DESeq2-normalised counts. Columns are sample library names; Cb = cerebellum, Hc = hippocampus.")),
  c(30, 120))

all_proc <- c(proc_files, quant_files)
md5 <- data.frame(`file name` = basename(all_proc), `file type` = ifelse(grepl("quant", all_proc), "Salmon quant", "matrix"),
                  `file checksum (MD5)` = unname(tools::md5sum(all_proc)), check.names = FALSE)
md5 <- rbind(md5, data.frame(`file name` = todo(paste0(rep(samples$np, each = 2), c(" R1", " R2"), " FASTQ")),
                             `file type` = "fastq", `file checksum (MD5)` = todo("md5sum of the raw file"),
                             check.names = FALSE))
add_sheet("MD5", md5, c(40, 14, 40))

saveWorkbook(wb, file.path(out_dir, "GEO_metadata_draft.xlsx"), overwrite = TRUE)
cat("Saved:", out_dir, "\n")
print(md5[seq_along(proc_files), ], row.names = FALSE)
