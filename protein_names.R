library(biomaRt)
library(KEGGREST)
library(dplyr)
library(rentrez)

# Lecture du fichier
data <- read.delim("all_counts_with_header.txt")

get_protein_name <- function(gene) {
  Sys.sleep(0.34)
  # Recherche du GeneID NCBI
  search_res <- entrez_search(db = "gene",
                              term = paste0(gene, " [Gene Name] AND Staphylococcus aureus [Organism]"))
  if(length(search_res$ids) == 0) return(NA)
  
  gene_id <- search_res$ids[1]
  
  # Récupérer les summaries pour trouver le RefSeq protéine
  summary <- entrez_summary(db="gene", id=gene_id)
  
  # Extraire le RefSeq protéine principal
  prot_ref <- summary$otheraliases  # parfois la liste des protéines (ou voir summary$genomicinfo)
  
  # Option plus fiable : on peut chercher le gène dans la base protein par [gene] AND [organism]
  prot_search <- entrez_search(db="protein",
                               term=paste0(gene, " [Gene] AND Staphylococcus aureus [Organism]"),
                               retmax=1)
  if(length(prot_search$ids)==0) return(NA)
  
  prot_summary <- entrez_summary(db="protein", id=prot_search$ids[1])
  prot_name = prot_summary$title
  if(!is.na(prot_name)) {
    # Supprimer le préfixe principal
    prot_name <- sub("^RecName: Full=", "", prot_name)
    # Remplacer les AltName par "/"
    prot_name <- gsub("; AltName: Full=", " / ", prot_name)
  }
  
  return(prot_name)  # contient le nom de la protéine
}

gene_info <- t(sapply(data$Geneid, get_protein_name))

# 1. Transformer la matrice en vecteur nommé
gene_info_vec <- as.vector(gene_info)
names(gene_info_vec) <- colnames(gene_info)

# 2. Créer un data frame avec ces infos
gene_info_df <- data.frame(
  Geneid = names(gene_info_vec),
  Protein_name = gene_info_vec,
  stringsAsFactors = FALSE
)

# 3. Joindre proprement à ton tableau de comptage
data_final <- merge(data, gene_info_df, by = "Geneid", all.x = TRUE)



