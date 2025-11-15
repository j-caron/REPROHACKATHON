data <- read.delim("all_counts_with_header.txt", row.names=1)
data_article = read.delim("GSE139659_IPvsctrl.complete.xls.gz")

# Suppose que ton dataframe s'appelle df
names(data)[1:3] <- paste0("per_", 1:3)
names(data)[4:6] <- paste0("control_", 1:3)

library(edgeR) #bioconductor, taper edgeR !
data_cpm = cpm(data)

library(DESeq2) #bioconductor aussi

cond = factor(c(rep("per",3),rep("control",3)))
cnts = data[rowSums(data)>10,]

dds = DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~cond)

dds = DESeq(dds)

res = results(dds)

res[order(res$pvalue)[1:100],] 


table(res$padj<.05) 
table(res$pvalue<.05/dim(cnts)[1])

par(scipen = 100)

hist(res$pvalue) 
par(bg = "gray80")
plotMA(res, ylim = c(-5, 5), colSig = 'red', colNonSig = 'black',xaxt = "n")
abline(h = pretty(res$log2FoldChange), col = "white", lty = 3)

seq = c(0,2,4,6)
powers <- 10^seq

abline(v = powers, col = "white", lty = 3)

axis(1, at = powers, labels = parse(text = paste0("10^",seq)))

res = as.data.frame(res)

data$Name = rownames(data)
res$Name = rownames(res)

data_perso = merge(data, res, by = "Name")

library(ggplot2)
library(ggrepel)

# -------------------
# Volcano plot
# -------------------
res$negLog10P <- -log10(res$pvalue)
res$significant <- "Non significatif"
res$significant[!is.na(res$pvalue) & res$pvalue < 0.05 & abs(res$log2FoldChange) >= 1] <- "Significatif"

volcano <- ggplot(res, aes(x = log2FoldChange, y = negLog10P, color = significant, label = gene)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("gray60", "red")) +
  labs(title = "Volcano plot - DESeq2",
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme_minimal(base_size = 14) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

# -------------------
# PCA
# -------------------
# Transformation variance stabilizing
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = "cond", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

center_x <- mean(pcaData$PC1)
center_y <- mean(pcaData$PC2)

pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = cond)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "PCA - DESeq2") +
  geom_segment(aes(x = center_x, y = center_y, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "gray50") +
  geom_text_repel(aes(label = rownames(pcaData)), max.overlaps = 20) +
  theme_minimal(base_size = 14)


ggsave("volcano_plot.png", volcano, width = 8, height = 6)
ggsave("PCA_plot.png", pca_plot, width = 8, height = 6)

comp = merge(data_perso,data_article, by = "Name")

library(dplyr)

colonnes_a_verifier <- c(
  "per_1", "per_2", "per_3",
  "control_1", "control_2", "control_3",
  "ctrl4", "ctrl5", "ctrl6",
  "IP1", "IP2", "IP3",
  "log2FoldChange.y", "pvalue.y"
)

# Supprimer lignes avec 0 dans au moins une colonne
comp <- comp[rowSums(comp[, colonnes_a_verifier] == 0, na.rm = FALSE) == 0, ]

# Supprimer les lignes avec NA dans ces colonnes
comp <- na.omit(comp)


comp$diff_pers1 = abs((comp$per_1-comp$IP1)/comp$IP1)*100
comp$diff_pers2 = abs((comp$per_2-comp$IP2)/comp$IP2)*100
comp$diff_pers3 = abs((comp$per_3-comp$IP3)/comp$IP3)*100
comp$diff_cont1 = abs((comp$control_1-comp$ctrl4)/comp$ctrl4)*100
comp$diff_cont2 = abs((comp$control_2-comp$ctrl5)/comp$ctrl5)*100
comp$diff_cont3 = abs((comp$control_3-comp$ctrl6)/comp$ctrl6)*100
comp$diff_logFC = abs((comp$log2FoldChange.x-comp$log2FoldChange.y)/comp$log2FoldChange.y)*100

colonnes_a_moyenner <- c("diff_pers1","diff_pers2","diff_pers3","diff_cont1",
                         "diff_cont2","diff_cont3","diff_logFC")
moyennes_colonnes <- colMeans(comp[, colonnes_a_moyenner], na.rm = TRUE)

comp_means <- as.data.frame(t(moyennes_colonnes))
