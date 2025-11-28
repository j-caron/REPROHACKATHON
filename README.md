# ReproHackathon – Reproduction des données de `mettre article`

## Présentation du projet

*Description*

------------------------------------------------------------------------

## Structure du dépôt - *Refaire en mieux quand on aura les infos !*

reprohackathon-project/\

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

## Installation et exécution

### 1. Clonage du dépôt

```bash
git clone https://github.com/tonpseudo/reprohackathon-project.git
cd reprohackathon-project
```

### 2. Lancement de l'image Docker

```bash
code pour lancer l'image
```

### 3. Exécution du workflow Snakemake 

```bash
code pour lancer l'exécution
```

Les résultats seront automatiquement générés dans le dossier `output/`.

## Description du workflow

Les workflows Snakemake sont à exécuter successivement :

1. **Téléchargement des fastq**
    Snakefile_download

2.  **Pipeline NGS**
    Snakefile_analysis

3.  **Analyse statistique et visualisations**
    Snakefile_R

## Résultats

Les résultats de l'analyse sont disponibles dans le dossier `output/`.

| Type de résultat | Fichier                               | Description                               |
|------------------|---------------------------------------|-------------------------------------------|
| Graphique        |                                       | MA plot des données totales               |
| Graphique        |                                       | MA plot des gènes liés à la traduction    |
| Graphique        |                                       | Volcano plot                              |
| Graphique        |                                       | PCA plot                                  |

## Auteurs

**Là on met nos noms et ce qu'on a fait**
