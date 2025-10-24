
# Fichier de configuration
configfile: "config.yaml"

# Liste des identifiants SRA
SAMPLES = config["sra_accessions"]

rule all:
    input:
        expand("{mapdir}/{sample}.bam",
               mapdir=config["dataMapped"],
               sample=config["sra_accessions"])

# Téléchargement des fichiers fastq à partir des identifiant SRA
rule download_fastq:
    output:
        "{outdir}/{sample}.fastq"
    params:
        outdir = config["dataFastq"]
    container:
        "quay.io/biocontainers/sra-tools:3.0.0--pl5321h9f0ad1d_1"
    threads:
        config["threads"]
    shell:
        """
        mkdir -p {params.outdir}
        fasterq-dump {wildcards.sample} -O {params.outdir} --threads {threads}
        """

# Téléchargement du génome 
rule download_genome:
    output:
        fasta = "{genomedir}/genome.fasta"
    params:
        genomedir = config["dataGenome"],
        genome_url = config["genomeURL"]
    container:
        "quay.io/biocontainers/wget:1.21.3--h0b77cf5_0"
    shell:
        """
        mkdir -p {params.genomedir}
        wget -O {params.genomedir}/genome.fna.gz {params.genome_url}
        gunzip -c {params.genomedir}/genome.fna.gz > {output.fasta}
        rm {params.genomedir}/genome.fna.gz
        """

rule trimming:
    input:
        "{fastqdir}/{sample}.fastq"
    output:
        "{trimdir}/{sample}_trimmed.fq.gz"
    params:
        fastqdir = config["dataFastq"],
        trimdir = config["dataTrimmed"]
    container:
        "quay.io/biocontainers/cutadapt:1.11--py35_0"
    threads:
        config["threads"]
    shell:
        # Même paramètre que l'article 
        """
        mkdir -p {params.trimdir}
        cutadapt -q 20 -m 25 -o {output} {input}   
        """

rule index_genome:
    input:
        fasta = rules.download_genome.output.fasta
    output:
        expand("{genomedir}/genome.{ext}",
               genomedir=config["dataGenome"],
               ext=["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"])
    container:
        # Dificulté à créer une image avec la même version bowtie que l'article
        # Utilisation de la version la plus proche bowtie1.0.0
        "quay.io/biocontainers/bowtie:1.0.0--py27_0"
    shell:
        # Même paramètre que l'article (par défault)
        """
        bowtie-build {input.fasta} {config[dataGenome]}/genome
        """

rule mapping:
    input:
        trimmed = "{trimdir}/{sample}_trimmed.fq.gz",
        index = expand("{genomedir}/genome.{ext}",
                       genomedir=config["dataGenome"],
                       ext=["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"])
    output:
        "{mapdir}/{sample}.bam"
    params:
        mapdir = config["dataMapped"]
    container:
        "quay.io/biocontainers/bowtie:1.0.0--py27_0"
    threads:
        config["threads"]
    shell:
        """
        mkdir -p {params.mapdir}
        bowtie -x {config[dataGenome]}/genome -U {input.trimmed} -p {threads} \
        | samtools view -Sb - > {output}
        """

