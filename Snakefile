
# Fichier de configuration
configfile: "config.yaml"

# Liste des identifiants SRA
SAMPLES = config["sra_accessions"]

rule all:
    input:
        expand("{outdir}/{sample}.fastq",
               outdir=config["dataFastq"],
               sample=SAMPLES)

# Téléchargement des fichiers fastq à partir des identifiant SRA
rule download:
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
