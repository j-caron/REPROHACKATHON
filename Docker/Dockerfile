# TOOLS NEEDED
# sra-toolkit:latest 
# fastqc:latest
# cutadapt:1.11
# TrimGalore:latest
# Bowtie:0.12.7
# samtools:latest
# Feature count:1.4.6-p3

FROM ubuntu:18.04 
# Ubuntu 18.04 is used because newer versions no longer include Python 2.7, required for Cutadapt 1.11.

LABEL description="ToolsForReprohackatonProject"

RUN apt-get update && apt-get install -y wget unzip tar build-essential python python-dev python-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /opt 

# Cutadapt
RUN pip install cutadapt==1.11 

# For the remaining tools, the archives were downloaded from official sources, extracted, 
# and all executables were moved to /usr/local/bin/ to make them accessible system-wide.

# SRA Toolkit
RUN wget -t 5 --wait=15 --retry-connrefused -O sra-toolkit.tar.gz "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-ubuntu64.tar.gz" && \
    tar -xzf sra-toolkit.tar.gz && \
    rm sra-toolkit.tar.gz && \
    chmod +x sratoolkit.3.2.1-ubuntu64/bin/* && \
    mv sratoolkit.3.2.1-ubuntu64/bin/* /usr/local/bin/

# Add SRA Toolkit to PATH
ENV PATH="/opt/sratoolkit.3.2.1-ubuntu64/bin:${PATH}"

# Java (required for FastQC)
RUN apt-get update && apt-get install -y default-jre && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# FastQC 
RUN wget -q -O fastqc.zip "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip" && \
    unzip fastqc.zip && \
    rm fastqc.zip && \
    chmod +x FastQC/fastqc && \
    mv FastQC /usr/local/FastQC && \
    ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc

# Trim Galore 
RUN wget -q -O TrimGalore.tar.gz "https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.10.tar.gz" && \
    tar -xzf TrimGalore.tar.gz && \
    rm TrimGalore.tar.gz && \
    chmod +x TrimGalore-0.6.10/trim_galore && \
    mv TrimGalore-0.6.10/trim_galore /usr/local/bin/

# Bowtie
RUN wget -q -O bowtie.zip "https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip/download" && \
    unzip bowtie.zip && \
    rm bowtie.zip && \
    chmod +x bowtie-0.12.7/bowtie* && \
    mv bowtie-0.12.7/bowtie* /usr/local/bin/

# Samtools (compiled from source)
RUN apt-get update && apt-get install -y \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    zlib1g-dev && \
    wget -q -O samtools.tar.bz2 "https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2" && \
    tar -xjf samtools.tar.bz2 && rm samtools.tar.bz2 && \
    cd samtools-1.22.1 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && make install && \
    cd .. && rm -rf samtools-1.22.1

# Subread (includes FeatureCounts)
RUN wget -q -O subread.tar.gz "https://sourceforge.net/projects/subread/files/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz/download" && \
    tar -xzf subread.tar.gz && \
    rm subread.tar.gz && \
    chmod -R +x subread-1.4.6-p3-Linux-x86_64/bin && \
    mv subread-1.4.6-p3-Linux-x86_64/bin/* /usr/local/bin/


ENV PATH="/usr/local/bin:${PATH}"


WORKDIR /data
