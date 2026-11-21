library(DESeq2)
library(ggplot2)

PCA_chart_drawer <- function(vsd, groups, tissue, group_colors, output_file) {
    
  pca <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
  pct <- round(100 * attr(pca, "percentVar"))
  pca$label <- gsub("_", " ", pca$name)
  p_pca <- ggplot(pca, aes(PC1, PC2, color = group, label = label)) +
    geom_point(size = 4) +
    geom_text(size = 3, vjust = -1.2, show.legend = FALSE) +
    scale_color_manual(values = group_colors, breaks = groups, labels = gsub("_", " ", groups)) +
    xlab(paste0("PC1: ", pct[1], "% variance")) +
    ylab(paste0("PC2: ", pct[2], "% variance")) +
    labs(title = paste0("PCA: ", paste(groups, collapse = " vs "), " (", tissue, ")"), color = "Group") +
    theme_bw()
  ggsave(output_file, p_pca, width = 7.5, height = 5)
}