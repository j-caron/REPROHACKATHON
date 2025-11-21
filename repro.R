data = read.table("GSE139659_IPvsctrl.complete.txt",
                   header = TRUE,      # première ligne = noms de colonnes
                   sep = "\t",         # séparateur = tabulation
                   check.names = FALSE) # garde les noms tels quels
data1 =  read.table("all_counts.txt", header = TRUE,      # première ligne = noms de colonnes
                    sep = "\t",         # séparateur = tabulation
                    check.names = FALSE)

data2  =  read.table("all_counts_with_header.txt", header = TRUE,      # première ligne = noms de colonnes
                            sep = "\t",         # séparateur = tabulation
                            check.names = FALSE)
data_projet = data[,1:11]

data2 = read.delim("~/GSE139659_IPvsctrl.complete.txt", row.names=2)
data2$significant = data2$padj < 0.05
data2 <- read.table("GSE139659_IPvsctrl.complete.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data = read.table("all_counts_with_header.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data3 = read.table("all_counts.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
table <- read.table("kegg_pathways.tsv", sep = "\t", header = FALSE)

valeurs <- c("path:sao03012", "path:sao03029", "path:sao03010", "path:sao00970")

table_filtre <- table[ table$V2 %in% valeurs , ]
table$V1 <- sub("^sao:", "", table$V1)
data2$gene_id = rownames(data2)
match_result <- data2$gene_id %in% table_filtre$V1
prop.table(table(match_result))
datatrans = data2[match_result,]
datatrans_clean <- datatrans[complete.cases(datatrans[, c("baseMean", "log2FoldChange", "significant")]), ]

library(ggplot2)

ggplot(datatrans_clean, aes(x = log2(baseMean), y = log2FoldChange,
                            color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  
  # Axe X : 0 → 20 par pas de 2
  scale_x_continuous(limits = c(0, 20),
                     breaks = seq(0, 20, 2)) +
  
  # Axe Y : -6 → 5 par pas de 1
  scale_y_continuous(limits = c(-6, 5),
                     breaks = seq(-6, 5, 1)) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Log2 base Mean",
       y = "Log2 Fold Change",
       color = "Significant")


