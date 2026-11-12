# ============================================
# Exercise_GF: barplot del numero di geni differenzialmente espressi (su/giu')
# per i tre confronti del modello a 3 gruppi, cervelletto e ippocampo.
# Legge le tabelle create da three_groups_exercise_gf.R (non rilancia DESeq2).
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_barplot.R")
# ============================================
library(ggplot2)

deseq_dir <- "results/rna-seq/Exercise_GF/deseq2"
padj_cut  <- 0.05
tissues   <- c(Cb = "Cerebellum", Hc = "Hippocampus")

contrasts <- data.frame(
  label = c("Bex3KO vs WT Control", "Bex3KO EE vs WT Control", "Bex3KO EE vs Bex3KO"),
  lfc   = c("log2FC_Bex3KO_vs_WT", "log2FC_Bex3KO_EE_vs_WT", "log2FC_Bex3KO_EE_vs_Bex3KO"),
  padj  = c("padj_Bex3KO_vs_WT", "padj_Bex3KO_EE_vs_WT", "padj_Bex3KO_EE_vs_Bex3KO")
)

# Conteggio dei geni su e giu' per ogni tessuto e confronto
counts <- do.call(rbind, lapply(names(tissues), function(t) {
  tab <- read.csv(file.path(deseq_dir, paste0("three_groups_", t), "three_groups_results.csv"))
  do.call(rbind, lapply(seq_len(nrow(contrasts)), function(i) {
    sig <- !is.na(tab[[contrasts$padj[i]]]) & tab[[contrasts$padj[i]]] < padj_cut
    lfc <- tab[[contrasts$lfc[i]]]
    data.frame(tissue = tissues[[t]], contrast = contrasts$label[i],
               direction = c("Up", "Down"),
               n = c(sum(sig & lfc > 0), sum(sig & lfc < 0)))
  }))
}))
counts$contrast  <- factor(counts$contrast, levels = contrasts$label)
counts$direction <- factor(counts$direction, levels = c("Up", "Down"))
counts$tissue    <- factor(counts$tissue, levels = tissues)
# Barre "giu'" disegnate verso il basso
counts$y <- ifelse(counts$direction == "Up", counts$n, -counts$n)

deg_barplot <- function(d, title) {
  ggplot(d, aes(contrast, y, fill = direction)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = n, vjust = ifelse(direction == "Up", -0.4, 1.4)), size = 3.5) +
    geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
    scale_y_continuous(labels = abs, expand = expansion(mult = 0.12)) +
    labs(title = title,
         subtitle = paste0("Three-group DESeq2 model, padj < ", padj_cut,
                           "\nUp/Down = higher/lower in the first group of each comparison"),
         x = NULL, y = "Number of DEGs", fill = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

# Un grafico (e la tabella dei numeri) per tessuto, nella sua cartella
for (t in names(tissues)) {
  d <- counts[counts$tissue == tissues[[t]], ]
  out_dir <- file.path(deseq_dir, paste0("three_groups_", t))
  p <- deg_barplot(d, paste0("Differentially expressed genes: ", tissues[[t]]))
  ggsave(file.path(out_dir, "deg_barplot.png"), p, width = 7, height = 5.5)
  write.csv(d[, c("tissue", "contrast", "direction", "n")],
            file.path(out_dir, "deg_counts.csv"), row.names = FALSE)
}

print(counts[, c("tissue", "contrast", "direction", "n")], row.names = FALSE)
