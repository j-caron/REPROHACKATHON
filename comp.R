#Script pour la comparaison des tableaux de données
#input = le tableau final comp de la fin de DESeq2.R

# => à voir si c'est intéressant, bien, tout ça, vous me dites

library(dplyr)

colonnes_a_verifier <- c(
  "per_1", "per_2", "per_3",
  "control_1", "control_2", "control_3",
  "ctrl4", "ctrl5", "ctrl6",
  "IP1", "IP2", "IP3"
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

colonnes_a_moyenner <- c("diff_pers1","diff_pers2","diff_pers3","diff_cont1",
                         "diff_cont2","diff_cont3","diff_logFC")
moyennes_colonnes <- colMeans(comp[, colonnes_a_moyenner], na.rm = TRUE)

comp_means <- as.data.frame(t(moyennes_colonnes))