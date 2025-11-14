data <- read.delim("~/all_counts_with_header.txt", row.names=1)

head(data)
summary(data)

boxplot(data)

# Suppose que ton dataframe s'appelle df
names(data)[1:3] <- paste0("per_", 1:3)
names(data)[4:6] <- paste0("control_", 1:3)



library(edgeR) #bioconductor, taper edgeR !
data_cpm = cpm(data)

library(DESeq2) #bioconductor aussi

####Peut-on utiliser une loi de Poisson pour analyser l'expression différentielle ?####
vectMean = apply(data,1,mean)
vectVar = apply(data,1,var)
plot(log(vectMean+1),log(vectVar+1))
abline(0,1)

cond = factor(c(rep("per",3),rep("control",3)))
cnts = data[rowSums(data)>10,]

dds = DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~cond)

dds = DESeq(dds)

res = results(dds)
#LogFC et p-values

res[order(res$pvalue)[1:100],] #Plus petits p-values


table(res$padj<.05) 
table(res$pvalue<.05/dim(cnts)[1])

par(scipen = 100)

hist(res$pvalue) 
par(bg = "gray80")
plotMA(res, ylim = c(-5, 5), colSig = 'red', colNonSig = 'black',xaxt = "n")
abline(h = pretty(res$log2FoldChange), col = "white", lty = 3)

# puissances de 10 adaptées à tes données
seq = c(0,2,4,6)
powers <- 10^seq   # 1, 10, 100, 1000, 10000, 100000, 1e6

abline(v = powers, col = "white", lty = 3)

axis(1, at = powers, labels = parse(text = paste0("10^",seq)))




res = as.data.frame(res)

data$gene = rownames(data)
res$gene = rownames(res)

data_perso = merge(data, res, by = "gene")

data_article = read.delim("~/GSE139659_IPvsctrl.complete.xls.gz")

