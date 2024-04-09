#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
NAMES = ["PacBio.full", "Ultralong.full", "NSS.full"]
CALLERS = ["clair2", "longshot", "nanocaller"][2:]

rule all:
    input:
        expand(OUTDIR + "/stratification/{caller}/{name}", caller=CALLERS, name=NAMES),

rule happy:
    input:
        vcf1 = SNP_VCF,
        vcf2 = OUTDIR + "/snvs/concated/{caller}/{name}.vcf.gz",
        tsv = STRATIFICATION_TSV,
        fasta = FASTA
    output:
        directory(OUTDIR + "/stratification/{caller}/{name}")
    log:
        OUTDIR + "/stratification/{caller}/{name}.log"
    conda:
        "happy"
    threads:
        8
    shell:
        """(
        mkdir -p {output}
        hap.py -r {input.fasta} -o {output}/benchmark --threads {threads} \
            --engine xcmp --no-roc --stratification {input.tsv} \
            {input.vcf1} {input.vcf2} ) &> {log}
        """
