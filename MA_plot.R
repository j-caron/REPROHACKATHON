#!/usr/bin/env Rscript
library(methods)
library(ggplot2)
options(bitmapType = "cairo")


data_perso_file <- snakemake@input[["data_perso"]]
kegg_file       <- snakemake@input[["kegg"]]
out_file        <- snakemake@output[[1]]

# -----------------------------
# Charger data_perso.csv
# -----------------------------
data2 <- read.csv(data_perso_file, stringsAsFactors = FALSE)

# Recréer colonne significant (n'existait pas dans data_perso)
data2$significant <- "Non significatif"
data2$significant[!is.na(data2$padj) &
                  data2$padj < 0.05 &
                  abs(data2$log2FoldChange) >= 1] <- "Significatif"

# Ajouter gene_id
data2$gene_id <- data2$Name

# -----------------------------
# Charger KEGG
# -----------------------------
table <- read.table(kegg_file, sep = "\t", header = FALSE)

# Pathways ciblés
valeurs <- c("path:sao03012", "path:sao03029", "path:sao03010", "path:sao00970")

table_filtre <- table[table$V2 %in% valeurs, ]
table_filtre$V1 <- sub("^sao:", "", table_filtre$V1)

# -----------------------------
# Filtrer data_perso avec KEGG
# -----------------------------
match_result <- data2$gene_id %in% table_filtre$V1
datatrans <- data2[match_result, ]

# Vérifie colonnes
stopifnot(all(c("baseMean","log2FoldChange","significant") %in% colnames(datatrans)))

datatrans_clean <- datatrans[
    complete.cases(datatrans[, c("baseMean","log2FoldChange","significant")]),
]

# -----------------------------
# MA Plot
# -----------------------------
p <- ggplot(datatrans_clean, aes(x = log2(baseMean),
                                 y = log2FoldChange,
                                 color = significant)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("grey", "red")) +
    scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 2)) +
    scale_y_continuous(limits = c(-6, 5), breaks = seq(-6, 5, 1)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(x = "Log2 baseMean", y = "Log2 Fold Change")

ggsave(out_file, p, width = 8, height = 6)
