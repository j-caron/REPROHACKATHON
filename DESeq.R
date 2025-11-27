#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Snakemake interface : récupère les fichiers du Snakefile
# ------------------------------------------------------------------------------
counts_file  <- snakemake@input[["counts"]]
article_file <- snakemake@input[["article"]]
output_dir   <- snakemake@params[["outdir"]]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Chargement des données
# ------------------------------------------------------------------------------
data <- read.delim(counts_file, row.names = 1)
data_article <- read.delim(article_file)

names(data)[1:3] <- paste0("per_", 1:3)
names(data)[4:6] <- paste0("control_", 1:3)

library(edgeR)
data_cpm = cpm(data)

library(DESeq2)

cond = factor(c(rep("per",3),rep("control",3)))
cnts = data[rowSums(data)>10,]

dds = DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~cond)
dds = DESeq(dds)
res = results(dds)
res = as.data.frame(res)

# ------------------------------------------------------------------------------
# Fusion
# ------------------------------------------------------------------------------
data$Name = rownames(data)
res$Name = rownames(res)
data_perso = merge(data, res, by = "Name")

library(ggplot2)
library(ggrepel)
options(bitmapType = "cairo")


# Volcano
res$negLog10P <- -log10(res$pvalue)
res$significant <- "Non significatif"
res$significant[!is.na(res$pvalue) & res$pvalue < 0.05 & abs(res$log2FoldChange) >= 1] <- "Significatif"

volcano <- ggplot(res, aes(x = log2FoldChange, y = negLog10P, color = significant)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal()

# PCA
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "cond", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = cond)) +
  geom_point(size=4) +
  theme_minimal()

# ------------------------------------------------------------------------------
# Sauvegardes
# ------------------------------------------------------------------------------
ggsave(file.path(output_dir, "volcano_plot.png"), volcano, width = 8, height = 6)
ggsave(file.path(output_dir, "PCA_plot.png"), pca_plot, width = 8, height = 6)
write.csv(data_perso, file.path(output_dir, "data_perso.csv"))

# Comparaison
comp = merge(data_perso,data_article, by="Name")
write.csv(comp, file.path(output_dir, "comparison_means.csv"))
