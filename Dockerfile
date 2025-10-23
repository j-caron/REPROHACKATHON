# Image de base : Debian stable
FROM ubuntu:22.04

LABEL maintainer="agash.uthayakumar0@gmail.com"
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /project

#Installation des dépendances système
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    curl \
    git \
    unzip \
    default-jre \
    ca-certificates \
    python3 \
    python3-pip \
    r-base \
    libbz2-dev liblzma-dev libncurses5-dev zlib1g-dev libssl-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Installation de Miniconda pour les outils bioinfo
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh \
    && /opt/conda/bin/conda init bash

#Création de l’environnement conda et installation des outils RNA-seq (alignement + comptage)
RUN /opt/conda/bin/conda config --system --remove channels defaults \
    && /opt/conda/bin/conda config --system --add channels conda-forge \
    && /opt/conda/bin/conda config --system --add channels bioconda \
    && /opt/conda/bin/conda create -y -n rnaseq_env python=3.10 mamba \
    && /opt/conda/bin/conda clean -afy


SHELL ["/bin/bash", "-lc"]

RUN source /opt/conda/bin/activate rnaseq_env \
    && mamba install -y -c bioconda -c conda-forge \
       hisat2 \
       sra-tools\
       bowtie2 \
       samtools \
       subread \ 
       bedtools \
       biopython \
       bwa \
       fastqc \
       multiqc \
       htseq\
    && mamba clean -afy

# Variables d’environnement
ENV CONDA_DEFAULT_ENV=rnaseq_env
ENV CONDA_PREFIX=/opt/conda/envs/rnaseq_env
ENV PATH=$CONDA_PREFIX/bin:$PATH

CMD ["/bin/bash"]
