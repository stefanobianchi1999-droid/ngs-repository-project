# ============================================
# Exercise_GF: PCA con ellissi attorno a ogni gruppo
#   - una PCA per tessuto (cervelletto, ippocampo) con i 3 gruppi
#   - una PCA di tutti i 18 campioni (gruppo x tessuto)
# Le ellissi RACCHIUDONO i campioni di ogni gruppo (come nei report del laboratorio):
# con 3 animali per gruppo un'ellisse di confidenza al 95% non si puo' calcolare.
# Stessa trasformazione dell'analisi DESeq2: VST (blind) e 500 geni piu' variabili.
#
# Da lanciare con la working directory nella root del progetto:
#   setwd("C:/Users/hp/Documents/GitHub/ngs-repository-project")
#   source("pipelines/rna-seq/scripts/r/three_groups_pca.R")
# ============================================
source("pipelines/rna-seq/scripts/r/exercise_gf_setup.R")
library(ggrepel)   # etichette che non si sovrappongono

group_levels <- c("WT_Control", "Bex3KO", "Bex3KO_EE")
group_colors <- c(WT_Control = "grey40", Bex3KO = "firebrick", Bex3KO_EE = "steelblue")
group_labels <- c(WT_Control = "WT Control", Bex3KO = "Bex3KO", Bex3KO_EE = "Bex3KO EE")

# Ellisse che racchiude i punti di un gruppo: centro = media, assi dalla covarianza,
# scalata perche' il punto piu' lontano stia sul bordo, piu' un margine del 25%
enclosing_ellipse <- function(x, y, pad = 1.25, n = 100) {
  pts <- cbind(x, y)
  ctr <- colMeans(pts)
  S <- cov(pts)
  S <- S + diag(1e-6 * max(diag(S), 1), 2)   # evita matrici singolari (punti allineati)
  d <- sqrt(max(mahalanobis(pts, ctr, S)))
  e <- eigen(S)
  t <- seq(0, 2 * pi, length.out = n)
  circle <- rbind(cos(t), sin(t))
  shape <- e$vectors %*% diag(sqrt(pmax(e$values, 0))) %*% circle * d * pad
  data.frame(x = ctr[1] + shape[1, ], y = ctr[2] + shape[2, ])
}

ellipses_for <- function(pca, by) {
  do.call(rbind, lapply(split(pca, pca[[by]]), function(d) {
    el <- enclosing_ellipse(d$PC1, d$PC2)
    el[[by]] <- d[[by]][1]
    el
  }))
}

# VST di un sottoinsieme di campioni, con lo stesso filtro dell'analisi DESeq2
vst_for <- function(cd) {
  txi <- subset_txi(txi_all, cd$sample)
  dds <- DESeqDataSetFromTximport(txi, colData = cd, design = ~ 1)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  dds <- dds[!chrom_of(rownames(dds)) %in% "Y", ]
  vst(dds, blind = TRUE)
}

# ---------- Una PCA per tessuto ----------
tissues <- c(Cb = "Cerebellum", Hc = "Hippocampus")
for (t in names(tissues)) {
  cd <- samples[samples$tissue == t, ]
  cd$group <- factor(cd$group, levels = group_levels)
  vsd <- vst_for(cd)
  pca <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
  pct <- round(100 * attr(pca, "percentVar"))
  pca$group <- factor(cd[pca$name, "group"], levels = group_levels)
  pca$label <- gsub("_", " ", pca$name)
  el <- ellipses_for(pca, "group")

  p <- ggplot(pca, aes(PC1, PC2, color = group)) +
    geom_polygon(data = el, aes(x, y, fill = group), alpha = 0.12, linewidth = 0.7) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = label), size = 3.4, box.padding = 0.6, min.segment.length = 0,
                    seed = 1, show.legend = FALSE) +
    scale_color_manual(values = group_colors, labels = group_labels) +
    scale_fill_manual(values = group_colors, labels = group_labels) +
    xlab(paste0("PC1: ", pct[1], "% variance")) +
    ylab(paste0("PC2: ", pct[2], "% variance")) +
    labs(title = paste0("PCA: ", tissues[[t]]),
         subtitle = "VST, 500 most variable genes; ellipses enclose the samples of each group",
         color = "Group", fill = "Group") +
    theme_bw() + theme(plot.title.position = "plot")
  ggsave(file.path(deseq_dir, paste0("three_groups_", t), "pca.png"), p, width = 8, height = 6)
  cat(tissues[[t]], ": PC1", pct[1], "% - PC2", pct[2], "%\n")
}

# ---------- PCA di tutti i 18 campioni ----------
cd <- samples
cd$group <- factor(cd$group, levels = group_levels)
vsd <- vst_for(cd)
pca <- plotPCA(vsd, intgroup = c("group", "tissue"), returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))
# plotPCA sovrascrive la colonna "group" con gruppo:tessuto -> la ricostruisco
pca$Group  <- factor(cd[pca$name, "group"], levels = group_levels)
pca$Tissue <- factor(tissues[cd[pca$name, "tissue"]], levels = tissues)
pca$cond   <- paste(pca$Group, pca$Tissue)
el <- ellipses_for(pca, "cond")
el$Group <- factor(sub(" .*", "", el$cond), levels = group_levels)

p <- ggplot(pca, aes(PC1, PC2, color = Group)) +
  geom_polygon(data = el, aes(x, y, group = cond, fill = Group), alpha = 0.10, linewidth = 0.6) +
  geom_point(aes(shape = Tissue), size = 4) +
  scale_color_manual(values = group_colors, labels = group_labels) +
  scale_fill_manual(values = group_colors, labels = group_labels) +
  xlab(paste0("PC1: ", pct[1], "% variance")) +
  ylab(paste0("PC2: ", pct[2], "% variance")) +
  labs(title = "PCA: all samples",
       subtitle = "VST, 500 most variable genes; ellipses enclose each group within each tissue",
       color = "Group", fill = "Group", shape = "Tissue") +
  theme_bw() + theme(plot.title.position = "plot")
ggsave(file.path(deseq_dir, "pca_all_samples.png"), p, width = 8, height = 6)
cat("All samples: PC1", pct[1], "% - PC2", pct[2], "%\n")
