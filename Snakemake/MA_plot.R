#!/usr/bin/env Rscript

# ============================================================================
# Snakemake inputs / outputs
# ============================================================================
data_file   <- snakemake@input[["data_perso"]]
kegg_file   <- snakemake@input[["kegg"]]
output_dir  <- snakemake@params[["outdir"]]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Chargement des librairies
# ============================================================================
library(ggplot2)
library(ggrepel)
library(KEGGREST)
library(stringr)

options(bitmapType = "native")

# ============================================================================
# Chargement des données
# ============================================================================
df <- read.csv(data_file, row.names = 1)
kegg_table <- read.table(kegg_file, sep = "\t", header = FALSE)

# ============================================================================
# Configuration
# ============================================================================
genes_interest <- c("infA", "infB", "infC", "tsf", "pth", "frr")
translation_pathways <- c("sao03010","sao00970","sao03012",
                          "sao03015","sao03029","sao03008")
translation_pathways_v2 <- paste0("path:", translation_pathways)

# ============================================================================
# Extraction des pathways (PATHWAY ou BRITE)
# ============================================================================
convert_kegg_to_gene_symbol <- function(gene_id) {
  out <- NA
  tryCatch({
    Sys.sleep(1)
    res <- KEGGREST::keggGet(paste0("sao:", gene_id))
    if (!is.null(res) && !is.null(res[[1]]$SYMBOL)) {
      out <- res[[1]]$SYMBOL
    }
  }, error = function(e) {})
  return(out)
}

convert_many_for <- function(ids) {
  sapply(ids, convert_kegg_to_gene_symbol)
}
df$gene_symbol <- convert_many_for(df$Name)
df$gene_symbol[df$Name == "SAOUHSC_00475"] <- "pth"


df$pathway <- NA
df$is_AAsynthesis <- FALSE

for (i in seq_len(nrow(df))) {

  gene <- tryCatch(
    keggGet(paste0("sao:", df$Name[i])),
    error = function(e) NULL
  )

  # ---- Cas 1 : PATHWAY trouvé ----
  if (!is.null(gene) && !is.null(gene[[1]]$PATHWAY)) {

    pw <- names(gene[[1]]$PATHWAY)
    df$pathway[i] <- paste(pw, collapse = ";")

    # détecter AA-tRNA biosynthesis
    if (length(pw) > 0 &&
        gene[[1]]$PATHWAY[[1]] == "Aminoacyl-tRNA biosynthesis") {
      df$is_AAsynthesis[i] = TRUE
    }

  } else {

    # ---- Cas 2 : récupérer via BRITE ----
    brite <- gene[[1]]$BRITE

    if (!is.null(brite)) {
      sao_matches <- regmatches(
        brite,
        gregexpr("BR:(sao[0-9]+)", brite)
      )

      sao_ids <- unique(unlist(lapply(sao_matches, function(x) {
        if (length(x) > 0) sub("BR:", "", x) else NULL
      })))

      if (length(sao_ids) > 0) {
        df$pathway[i] <- paste(sao_ids, collapse = ";")
      } else {
        df$pathway[i] <- NA
      }
    } else {
      df$pathway[i] <- NA
    }
  }
}

# ============================================================================
# Détection des gènes translationnels (KEGGREST + fichier KEGG)
# ============================================================================
df$is_translation <- sapply(df$pathway, function(x) {
  pw <- unlist(strsplit(as.character(x), ";"))
  any(pw %in% translation_pathways)
})

# fichier kegg_pathways.tsv
table_filtre <- kegg_table[kegg_table$V2 %in% translation_pathways_v2, ]
table_filtre$V1 <- sub("^sao:", "", table_filtre$V1)

genes_keggfile  <- unique(table_filtre$V1)
genes_keggrest  <- df[df$is_translation, ]$Name
genes_final     <- union(genes_keggrest, genes_keggfile)

df$is_translation_combined <- df$Name %in% genes_final

df_translation <- df[df$is_translation_combined, ]

# ============================================================================
# Détection AA-tRNA biosynthesis dans subset translation
# ============================================================================
df_translation$is_AAsynthesis <- FALSE

for (i in seq_len(nrow(df_translation))) {

  gene_data <- tryCatch(
    keggGet(paste0("sao:", df_translation$Name[i])),
    error = function(e) NULL
  )

  if (!is.null(gene_data)) {
    pw <- gene_data[[1]]$PATHWAY
    if (!is.null(pw) && length(pw) > 0) {
      if (pw[[1]] == "Aminoacyl-tRNA biosynthesis") {
        df_translation$is_AAsynthesis[i] = TRUE
      }
    }
  }
}

# ============================================================================
# Préparation du plot
# ============================================================================
df_translation$color_group <- ifelse(df_translation$padj < 0.05,
                                     "Significatif", "Non significatif")
df_translation$color_group <- factor(df_translation$color_group,
                                     levels = c("Non significatif", "Significatif"))

df_translation$shape_group <- ifelse(df_translation$is_AAsynthesis,
                                     "AA-tRNA synthetases", "Autres")
df_translation$shape_group <- factor(df_translation$shape_group,
                                     levels = c("Autres", "AA-tRNA synthetases"))

df_translation$label <- ifelse(df_translation$gene_symbol %in% genes_interest,
                               df_translation$gene_symbol, NA)

cols <- c("Non significatif" = "grey70", "Significatif" = "red")

# ============================================================================
# MA plot
# ============================================================================
p <- ggplot(df_translation, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(
    aes(fill = color_group),
    shape = 21,
    color = "black",
    size = 3,
    alpha = 0.9,
    stroke = 0.3
  ) +
  geom_point(
    data = subset(df_translation, is_AAsynthesis),
    shape = 21,
    fill = NA,
    color = "black",
    size = 2,
    stroke = 1.7,
    show.legend = FALSE
  ) +
  geom_point(
    data = data.frame(
      x = NA, y = NA, Cat = "AA-tRNA synthetases"
    ),
    aes(x = x, y = y, shape = Cat),
    size = 2,
    fill = "white",
    color = "black",
    stroke = 1.5
  ) +
  scale_fill_manual(values = cols, name = "Significativité") +
  scale_shape_manual(values = c("AA-tRNA synthetases" = 21),
                     name = "Catégorie") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text_repel(
    aes(label = label),
    size = 5,
    segment.size = 0.6,
    segment.color = "black",
    nudge_x = 0.6, nudge_y = 0.6
  ) +
  labs(x = "Log2 baseMean", y = "Log2 Fold Change") +
  theme_bw(base_size = 14)

# ============================================================================
# Sauvegarde
# ============================================================================
ggsave(file.path(output_dir, "MA_plot.png"), p, width = 10, height = 8)
