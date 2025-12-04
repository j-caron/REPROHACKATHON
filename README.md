# ReproHackathon – Reproduction des données d'un article scientifique

Reproduction des données de l'article "Intracellular Staphylococcus aureus persisters upon antibiotic exposure" *(Peyrusson et al., https://doi.org/10.1038/s41467-020-15966-7)*

## Structure du dépôt

. \
├── Docker \
│   ├── DockerFile_R \
│   ├── Dockerfile \
│   └── requirements.txt \
├── Outputs \
│   ├── counts \
│   │   ├── SRR10379721_counts.txt \
│   │   ├── SRR10379721_counts.txt.summary \
│   │   ├── SRR10379722_counts.txt \
│   │   ├── SRR10379722_counts.txt.summary \
│   │   ├── SRR10379723_counts.txt \
│   │   ├── SRR10379723_counts.txt.summary \
│   │   ├── SRR10379724_counts.txt \
│   │   ├── SRR10379724_counts.txt.summary \
│   │   ├── SRR10379725_counts.txt \
│   │   ├── SRR10379725_counts.txt.summary \
│   │   ├── SRR10379726_counts.txt \
│   │   ├── SRR10379726_counts.txt.summary \
│   │   ├── all_counts.txt \
│   │   └── all_counts_with_header.txt \
│   ├── fastq \
│   │   ├── SRR10379721.fastq \
│   │   ├── SRR10379722.fastq \
│   │   ├── SRR10379723.fastq \
│   │   ├── SRR10379724.fastq \
│   │   ├── SRR10379725.fastq \
│   │   └── SRR10379726.fastq \
│   ├── genome \
│   │   ├── genome.1.ebwt \
│   │   ├── genome.2.ebwt \
│   │   ├── genome.3.ebwt \
│   │   ├── genome.4.ebwt \
│   │   ├── genome.fasta \
│   │   ├── genome.gff \
│   │   ├── genome.rev.1.ebwt \
│   │   ├── genome.rev.2.ebwt \
│   │   └── kegg_pathways.tsv \
│   ├── mapped \
│   │   ├── SRR10379721.bam \
│   │   ├── SRR10379722.bam \
│   │   ├── SRR10379723.bam \
│   │   ├── SRR10379724.bam \
│   │   ├── SRR10379725.bam \
│   │   └── SRR10379726.bam \
│   ├── results \
│   │   ├── MA_plot.png \
│   │   ├── MA_plot_tot.png \
│   │   ├── PCA_plot.png \
│   │   ├── volcano_plot.png \
│   │   ├── data_perso.csv \
│   │   └── comparison_means.csv \
│   └── trimmed \
│       ├── SRR10379722_trimmed.fq.gz \
│       ├── SRR10379723_trimmed.fq.gz \
│       ├── SRR10379724_trimmed.fq.gz \
│       ├── SRR10379725_trimmed.fq.gz \
│       └── SRR10379726_trimmed.fq.gz \
├── Snakemake \
│   ├── Snakefile_analysis \
│   ├── Snakefile_download \
│   ├── Snakefile_R \
│   └── config.yaml \
├── arborescence.txt \
└── environment.yaml \

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
snakemake -s Snakemake/Snakefile_R --cores 16 --use-singularity -p
```

Les résultats seront automatiquement générés dans le dossier `output/`.\
Les images Docker viennent des Dockerfiles disponibles dans `Docker`, et sont disponibles sur le Docker Hub avec les liens _<docker://mmousg/reprohackaton-tools:v1.1>_, et _<docker://agash00/r341_deseq2:1.16>_.

## Description du workflow

Les workflows Snakemake sont les suivants :

1. **Snakefile_download**
   
Récupération des fastq et du génome de référence.

2.  **Snakefile_analysis**
   
Trimming des fastq avec cutadapt, indexation du génome avec bowtie, mapping des fastq sur le génome de référence avec bowtie, récupération du fichier BAM avec samtools et comptage avec featureCounts
  
3.  **Snakefile_R**
   
Analyse statistique des comptages obtenus avec DESeq2, ajout des id des gènes avec l'API KEGG et réalisation des graphes (MA plot, ACP, Volcano plot et MA plot des gènes liés à la traduction).

## Résultats

Les résultats de l'analyse sont disponibles dans le dossier `Output/results`.

| Type de résultat | Fichier                               | Description                               |
|------------------|---------------------------------------|-------------------------------------------|
| Graphique        | MA_plot_tot.png                       | MA plot des données totales               |
| Graphique        | MA_plot.png                           | MA plot des gènes liés à la traduction    |
| Graphique        | volcano_plot.png                      | Volcano plot                              |
| Graphique        | PCA_plot.png                          | PCA plot                                  |
| Tableau          | data_perso.csv                        | Comptages, logFC et nom des gènes         |
| Tableau          | comparion_means.csv                   | Comparaison des données de l'article et de nos données |

## Auteurs

CARON Jérémy, GUEYE Mouhamadou Moustapha, UTHAYAKUMAR Agash
