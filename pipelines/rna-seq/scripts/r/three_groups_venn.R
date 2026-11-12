# ============================================
# Exercise_GF: diagramma di Venn dei tre confronti del modello a 3 gruppi
# (geni con padj < 0.05), cervelletto e ippocampo, con legenda.
# Legge le tabelle create da three_groups_exercise_gf.R (non rilancia DESeq2).
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_venn.R")
# ============================================
library(ggplot2)

deseq_dir <- "results/rna-seq/Exercise_GF/deseq2"
padj_cut  <- 0.05
tissues   <- c(Cb = "Cerebellum", Hc = "Hippocampus")

contrasts <- data.frame(
  label = c("Bex3KO vs WT Control", "Bex3KO EE vs WT Control", "Bex3KO EE vs Bex3KO"),
  padj  = c("padj_Bex3KO_vs_WT", "padj_Bex3KO_EE_vs_WT", "padj_Bex3KO_EE_vs_Bex3KO"),
  color = c("steelblue", "forestgreen", "firebrick")
)

# Venn a 3 insiemi: cerchi + conteggio dei geni in ognuna delle 7 regioni
venn3_plot <- function(sets, title) {
  a <- sets[[1]]; b <- sets[[2]]; c <- sets[[3]]
  counts <- c(
    length(setdiff(a, union(b, c))), length(setdiff(b, union(a, c))), length(setdiff(c, union(a, b))),
    length(setdiff(intersect(a, b), c)), length(setdiff(intersect(a, c), b)),
    length(setdiff(intersect(b, c), a)), length(Reduce(intersect, list(a, b, c)))
  )
  legend_labels <- paste0(names(sets), " (", sapply(sets, length), " genes)")
  centers <- data.frame(x = c(-0.6, 0.6, 0), y = c(0.35, 0.35, -0.65))
  circle <- do.call(rbind, lapply(seq_len(3), function(i) {
    t <- seq(0, 2 * pi, length.out = 200)
    data.frame(x = centers$x[i] + cos(t), y = centers$y[i] + sin(t), set = legend_labels[i])
  }))
  circle$set <- factor(circle$set, levels = legend_labels)
  labels <- data.frame(
    x = c(-1.1, 1.1, 0, 0, -0.65, 0.65, 0),
    y = c(0.75, 0.75, -1.25, 0.85, -0.35, -0.35, 0.05),
    n = counts
  )
  ggplot() +
    geom_polygon(data = circle, aes(x, y, fill = set, color = set), alpha = 0.25, linewidth = 0.8) +
    geom_text(data = labels, aes(x, y, label = n), size = 5) +
    scale_fill_manual(values = setNames(contrasts$color, legend_labels)) +
    scale_color_manual(values = setNames(contrasts$color, legend_labels)) +
    coord_equal() + theme_void() +
    theme(legend.position = "bottom", legend.direction = "vertical",
          legend.title = element_blank(), legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = title,
         subtitle = paste0("DEGs (padj < ", padj_cut, ") per comparison; numbers = genes in each region"))
}

for (t in names(tissues)) {
  tab <- read.csv(file.path(deseq_dir, paste0("three_groups_", t), "three_groups_results.csv"))
  sets <- lapply(contrasts$padj, function(col) tab$gene_id[!is.na(tab[[col]]) & tab[[col]] < padj_cut])
  names(sets) <- contrasts$label
  p <- venn3_plot(sets, paste0("Differentially expressed genes: ", tissues[[t]]))
  ggsave(file.path(deseq_dir, paste0("three_groups_", t), "venn_pairwise.png"), p,
         width = 6.5, height = 7.5, bg = "white")
}
