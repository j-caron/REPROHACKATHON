# ReproHackathon – Analyse reproductible de données biologiques

## Présentation du projet

*Définition du projet –*

------------------------------------------------------------------------

## Structure du dépôt - *Refaire en mieux !*

reprohackathon-project/\

│\

├── Data/ \# Données brutes d'entrée (.csv, .tsv, etc.)\

├── Snakemake/ \# Snakefiles et règles Snakemake\

│ ├── Snakefiles/ \# Fichiers principaux du workflow\

│ └── rules/ \# Règles modulaires (étapes spécifiques)\

├── Docker/ \# Fichiers de construction de l'image Docker\

│ ├── Dockerfiles/ \# Images Docker du projet

│ └── requirements.txt \# Dépendances (Python, R, etc.)\

├── scripts/ \# Scripts d'analyse exécutés par Snakemake\

├── config/ \# Fichiers de configuration (paramètres, chemins)\

├── Output/ \# Résultats finaux (graphes, tableaux, modèles)\

└── README.md \# Documentation principale du projet\

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

Le workflow Snakemake est structuré en plusieurs étapes successives :

1. **Préparation des données**\
    Nettoyage, harmonisation et vérification des fichiers d'entrée.

2.  **Analyse principale**\
    Exécution des scripts d'analyse (R ou Python) définis dans les règles Snakemake.

3.  **Visualisation**\
    Production de figures et graphiques synthétiques enregistrés dans `output/`.

## Résultats

Les résultats de l'analyse sont disponibles dans le dossier `output/`.

| Type de résultat | Fichier exemple                       | Description                               |
|------------------|---------------------------------------|-------------------------------------------|
| Graphique        | `output/temperature_distribution.png` | MA plot des données totales               |
| Graphique        | `output/species_frequency.png`        | MA plot des gènes liés à la traduction    |
| Graphique        | `output/model_summary.txt`            | Volcano plot                              |
| Graphique        |                                       | PCA plot                                  |

## Auteurs

**Là on met nos noms et ce qu'on a fait**
