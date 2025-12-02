# ReproHackathon – Reproduction des données d'un article scientifique

Reproduction des données de l'article "Intracellular Staphylococcus aureus persisters upon antibiotic exposure" *(Peyrusson et al., https://doi.org/10.1038/s41467-020-15966-7)*

## Structure du dépôt - *Refaire en mieux quand on aura les infos !*

REPROHACKATHON/\

│\

├── Data/

├── Snakemake/

│ ├── Snakefiles/

├── Docker/

│ ├── Dockerfiles/

│ └── requirements.txt

├── config/

├── Output/

└── README.md 

---
## Environnement de calcul utilisé

Les analyses ont été réalisées sur une machine virtuelle fournie par l’IFB-Biosphère.

### 1. Caractéristiques de la VM
- 16 cœurs CPU
- 64 Go de RAM
- 120 Go de stockage
- Système d’exploitation : Ubuntu 22.04 (Appliance IFB)

### 2. Identifiant de l’appliance
- Nom : Appliance Ubuntu 22.04  
- ID : `ubuntu-ifb`  
- Recipe Git : https://gitlab.in2p3.fr/ifb-biosphere/apps/ubuntu-ifb


## Installation et exécution


### 1. Clonage du dépôt

```bash
git clone https://github.com/j-caron/REPROHACKATHON.git
cd REPROHACKATHON/
```

### 2. Pré-installations: Docker et Snakemake

Pour reproduire exactement l’environnement utilisé dans la VM :

```bash
conda env create -f environment.yaml
conda activate snakemake_env
```

### 3. Exécution des workflows Snakemake 

```bash
snakemake -s Snakemake/Snakefile_download --use-singularity -j 16 -p
snakemake -s Snakemake/Snakefile_analysis --use-singularity -j 8 -p
```

Les résultats seront automatiquement générés dans le dossier `output/`.\
Les images Docker viennent de *mettre liens*.

## Description du workflow

Les workflows Snakemake sont à exécuter successivement :

1. **Téléchargement des fastq**\
    *Snakefile_download*

2.  **Pipeline NGS**\
    *Snakefile_analysis*

3.  **Analyse statistique et visualisations**\
    *Snakefile_R*

## Résultats

Les résultats de l'analyse sont disponibles dans le dossier `Output/`.

| Type de résultat | Fichier                               | Description                               |
|------------------|---------------------------------------|-------------------------------------------|
| Graphique        |                                       | MA plot des données totales               |
| Graphique        |                                       | MA plot des gènes liés à la traduction    |
| Graphique        |                                       | Volcano plot                              |
| Graphique        |                                       | PCA plot                                  |

## Auteurs

CARON Jérémy, GUEYE Mouhamadou Moustapha, UTHAYAKUMAR Agash
