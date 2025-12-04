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

library(DESeq2)

cond = factor(c(rep("per",3),rep("control",3)))

dds = DESeqDataSetFromMatrix(data, DataFrame(cond), ~cond)
dds = DESeq(dds)
res = results(dds)
res_def = as.data.frame(res)

# ------------------------------------------------------------------------------
# Fusion
# ------------------------------------------------------------------------------
data$Name = rownames(data)
res_def$Name = rownames(res)
data_perso = merge(data, res_def, by = "Name")

library(ggplot2)
library(ggrepel)
options(bitmapType = "cairo")

# MA plot en base R

# Colonne significatif / non significatif
res_def$significant <- ifelse(
  !is.na(res_def$padj) & res_def$padj < 0.05,
  "Significatif",
  "Non significatif"
)

res_def$significant <- factor(res_def$significant, 
                             levels = c("Non significatif", "Significatif"))

# Couleurs
cols <- c("Non significatif" = "black", "Significatif" = "red")

# MA Plot
ma_plot <- ggplot(res_def, aes(x = log2(baseMean + 1), y = log2FoldChange)) +
  geom_point(aes(color = significant), alpha = 0.7, size = 2) +
  scale_color_manual(values = cols, name = "Significativité") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Log2(baseMean + 1)",
    y = "Log2 Fold Change",
    title = "MA Plot (ggplot2)"
  ) +
  theme_bw(base_size = 14)

ggsave(file.path(output_dir, "MA_plot_tot.png"),
       ma_plot, width = 10, height = 7)


# Volcano
res_def$negLog10P <- -log10(res$pvalue)
res_def$significant <- "Non significatif"
res_def$significant[!is.na(res$pvalue) &
                    res$pvalue < 0.05 &
                    abs(res$log2FoldChange) >= 1] <- "Significatif"

cols <- c("Non significatif" = "grey70",
          "Significatif" = "red")

volcano <- ggplot(res_def,
                  aes(x = log2FoldChange,
                      y = negLog10P,
                      color = significant)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = cols) +   # <-- obligatoire
  theme_minimal()

ggsave(file.path(output_dir, "volcano_plot.png"),
       volcano, width = 8, height = 6)


# PCA
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "cond", returnData = TRUE)

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = cond)) +
  geom_point(size=4) +
  theme_minimal()

ggsave(file.path(output_dir, "PCA_plot.png"), pca_plot, width = 8, height = 6)

#percentVar <- round(100 * attr(pcaData, "percentVar"))


# ------------------------------------------------------------------------------
# Sauvegardes
# ------------------------------------------------------------------------------
write.csv(data_perso, file.path(output_dir, "data_perso.csv"))

# Comparaison
# à refaire !
comp = merge(data_perso,data_article, by="Name")
write.csv(comp, file.path(output_dir, "comparison_means.csv"))

