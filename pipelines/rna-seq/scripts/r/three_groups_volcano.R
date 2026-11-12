# ============================================
# Exercise_GF: volcano plot dei tre confronti del modello a 3 gruppi.
# Per ogni tessuto:
#   - un volcano grande per ogni confronto (volcano_<confronto>.png)
#   - i tre confronti affiancati (volcano_three_contrasts.png)
# I geni piu' significativi sono etichettati con ggrepel (etichette che non si sovrappongono).
# Legge le tabelle create da three_groups_exercise_gf.R (non rilancia DESeq2).
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_volcano.R")
# ============================================
library(ggplot2)
library(ggrepel)

deseq_dir <- "results/rna-seq/Exercise_GF/deseq2"
padj_cut  <- 0.05
tissues   <- c(Cb = "Cerebellum", Hc = "Hippocampus")
lfc_max   <- 6    # oltre +-6 i geni (quasi sempre poco espressi) vengono disegnati sul bordo
n_labels  <- 15   # geni etichettati per volcano

contrasts <- data.frame(
  label = c("Bex3KO vs WT Control", "Bex3KO EE vs WT Control", "Bex3KO EE vs Bex3KO"),
  file  = c("Bex3KO_vs_WT", "Bex3KO_EE_vs_WT", "Bex3KO_EE_vs_Bex3KO"),
  lfc   = c("log2FC_Bex3KO_vs_WT", "log2FC_Bex3KO_EE_vs_WT", "log2FC_Bex3KO_EE_vs_Bex3KO"),
  padj  = c("padj_Bex3KO_vs_WT", "padj_Bex3KO_EE_vs_WT", "padj_Bex3KO_EE_vs_Bex3KO")
)
volcano_colors <- c(Up = "firebrick", Down = "steelblue", "Not significant" = "grey75")

# Tabella per il volcano di un confronto
volcano_data <- function(tab, i) {
  d <- data.frame(gene_name = tab$gene_name, chrom = tab$chrom, baseMean = tab$baseMean,
                  log2FoldChange = tab[[contrasts$lfc[i]]], padj = tab[[contrasts$padj[i]]],
                  contrast = contrasts$label[i])
  d <- d[!is.na(d$padj) & !is.na(d$log2FoldChange), ]
  # padj = 0 darebbe -log10 = Inf: lo fisso a 1e-300
  d$padj <- pmax(d$padj, 1e-300)
  # Geni fuori scala: li porto sul bordo dell'asse x e li disegno come triangoli
  d$out_of_range <- abs(d$log2FoldChange) > lfc_max
  d$log2FoldChange <- pmax(pmin(d$log2FoldChange, lfc_max), -lfc_max)
  d$signif <- ifelse(d$padj < padj_cut & d$log2FoldChange > 0, "Up",
              ifelse(d$padj < padj_cut & d$log2FoldChange < 0, "Down", "Not significant"))
  d[order(d$padj), ]
}

# Geni da etichettare: i piu' significativi, ben espressi, non fuori scala,
# esclusi i geni predetti (Gm..., ...Rik) che non dicono nulla a chi legge
label_genes <- function(d, n) {
  ok <- d$signif != "Not significant" & d$baseMean >= 100 & !d$out_of_range &
        !grepl("^Gm[0-9]|Rik$|^ENSMUSG", d$gene_name)
  head(d[ok, ], n)
}

for (t in names(tissues)) {
  out_dir <- file.path(deseq_dir, paste0("three_groups_", t))
  tab <- read.csv(file.path(out_dir, "three_groups_results.csv"))
  all_volc <- list(); all_lab <- list()

  # ---------- Un volcano grande per confronto ----------
  for (i in seq_len(nrow(contrasts))) {
    d   <- volcano_data(tab, i)
    lab <- label_genes(d, n_labels)
    all_volc[[i]] <- d; all_lab[[i]] <- label_genes(d, 8)
    n_up <- sum(d$signif == "Up"); n_down <- sum(d$signif == "Down")

    p <- ggplot(d, aes(log2FoldChange, -log10(padj), color = signif)) +
      geom_point(aes(shape = out_of_range), size = 1.6, alpha = 0.6) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
      geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey50") +
      geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
      geom_text_repel(data = lab, aes(label = gene_name), size = 4.5, fontface = "bold",
                      box.padding = 0.5, point.padding = 0.3, min.segment.length = 0,
                      max.overlaps = Inf, seed = 1, show.legend = FALSE) +
      annotate("label", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 5,
               label = paste0("Up: ", n_up, "   Down: ", n_down)) +
      coord_cartesian(xlim = c(-lfc_max, lfc_max)) +
      scale_color_manual(values = volcano_colors, breaks = names(volcano_colors)) +
      labs(title = paste0(contrasts$label[i], " (", tissues[[t]], ")"),
           subtitle = paste0("Colored = padj < ", padj_cut, "; log2FC > 0 = higher in ",
                             sub(" vs .*", "", contrasts$label[i]),
                             "; top ", n_labels, " significant genes labelled"),
           x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
      theme_bw(base_size = 15) +
      theme(plot.title.position = "plot", legend.position = "bottom",
            plot.title = element_text(face = "bold"))
    ggsave(file.path(out_dir, paste0("volcano_", contrasts$file[i], ".png")), p,
           width = 10, height = 8.5)
  }

  # ---------- I tre confronti affiancati ----------
  volc <- do.call(rbind, all_volc)
  lab  <- do.call(rbind, all_lab)
  volc$contrast <- factor(volc$contrast, levels = contrasts$label)
  lab$contrast  <- factor(lab$contrast, levels = contrasts$label)
  n_lab <- do.call(rbind, lapply(split(volc, volc$contrast), function(d)
    data.frame(contrast = d$contrast[1],
               label = paste0("Up: ", sum(d$signif == "Up"), "   Down: ", sum(d$signif == "Down")))))

  p <- ggplot(volc, aes(log2FoldChange, -log10(padj), color = signif)) +
    geom_point(aes(shape = out_of_range), size = 1.1, alpha = 0.6) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey50") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    geom_text_repel(data = lab, aes(label = gene_name), size = 3.8, fontface = "bold",
                    box.padding = 0.4, min.segment.length = 0, max.overlaps = Inf,
                    seed = 1, show.legend = FALSE) +
    geom_text(data = n_lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = -0.05, vjust = 1.6, size = 4.2) +
    facet_wrap(~ contrast, nrow = 1) +
    coord_cartesian(xlim = c(-lfc_max, lfc_max)) +
    scale_color_manual(values = volcano_colors, breaks = names(volcano_colors)) +
    labs(title = paste0("Volcano plots, three-group model: ", tissues[[t]]),
         subtitle = paste0("Colored = padj < ", padj_cut,
                           "; log2FC > 0 = higher in the first group of each comparison; ",
                           "triangles = |log2FC| > ", lfc_max, " (drawn at the edge)"),
         x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) +
    theme_bw(base_size = 14) +
    theme(plot.title.position = "plot", legend.position = "bottom",
          strip.text = element_text(size = 13, face = "bold"))
  ggsave(file.path(out_dir, "volcano_three_contrasts.png"), p, width = 18, height = 7.5)
}
