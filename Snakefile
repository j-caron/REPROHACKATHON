# Fichier de configuration
configfile: "config.yaml"

shell.executable("/bin/sh")

SAMPLES = config["sra_accessions"]

# Règle finale
rule all:
    input:
        expand(
            f"{config['dataMapped']}/{{sample}}.bam",
            sample=SAMPLES
        )

# Téléchargement des fichiers FASTQ depuis SRA
rule download_fastq:
    output:
        f"{config['dataFastq']}/{{sample}}.fastq"
    container:
        "docker://ncbi/sra-tools:latest"
    threads:
        config["threads"]
    shell:
        """
        mkdir -p {config[dataFastq]}
        fasterq-dump {wildcards.sample} -O {config[dataFastq]} --threads {threads}
        """

# Téléchargement du génome de référence
rule download_genome:
    output:
        f"{config['dataGenome']}/genome.fasta"
    params:
        genome_url = config["genomeURL"]
    container:
        "docker://netdata/wget:latest"
    shell:
        """
        mkdir -p {config[dataGenome]}
        wget -O {config[dataGenome]}/genome.fna.gz {params.genome_url}
        gunzip -c {config[dataGenome]}/genome.fna.gz > {output}
        rm {config[dataGenome]}/genome.fna.gz
        """

# Trimming des lectures (Cutadapt)
rule trimming:
    input:
        f"{config['dataFastq']}/{{sample}}.fastq"
    output:
        f"{config['dataTrimmed']}/{{sample}}_trimmed.fq.gz"
    container:
        "X"
    threads:
        config["threads"]
    shell:
        """
        mkdir -p {config[dataTrimmed]}
        cutadapt -q 20 -m 25 -o {output} {input}
        """

# Indexation du génome (Bowtie)
rule index_genome:
    input:
        f"{config['dataGenome']}/genome.fasta"
    output:
        expand(
            "{genomedir}/genome.{ext}",
            genomedir=config["dataGenome"],
            ext=["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt",
                 "rev.1.ebwt", "rev.2.ebwt"]
        )
    params:
        genomedir = config["dataGenome"]
    container:
        "X"
    shell:
        """
        bowtie-build {input} {params.genomedir}/genome
        """

# Mapping des lectures (Bowtie)
rule mapping:
    input:
        trimmed = f"{config['dataTrimmed']}/{{sample}}_trimmed.fq.gz",
        index = expand(
            "{genomedir}/genome.{ext}",
            genomedir=config["dataGenome"],
            ext=["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt",
                 "rev.1.ebwt", "rev.2.ebwt"]
        )
    output:
        f"{config['dataMapped']}/{{sample}}.bam"
    params:
        genomedir = config["dataGenome"],
        mapdir = config["dataMapped"]
    container:
        "X"
    threads:
        config["threads"]
    shell:
        """
        mkdir -p {params.mapdir}
        bowtie -x {params.genomedir}/genome -U {input.trimmed} -p {threads} \
        | samtools view -Sb - > {output}
        """
