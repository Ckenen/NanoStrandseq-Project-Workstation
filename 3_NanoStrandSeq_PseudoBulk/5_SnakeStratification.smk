#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"

names = ["PacBio.full", "Ultralong.full", "NSS.full"]
callers = ["clair2", "longshot", "nanocaller"][2:]

rule all:
    input:
        expand(outdir + "/stratification/{caller}/{name}", caller=callers, name=names),

rule happy:
    input:
        vcf1 = BENCHMARK_VCF,
        vcf2 = outdir + "/snvs/concated/{caller}/{name}.vcf.gz",
        tsv = "data/GRCh38_stratification/v3.1-GRCh38-all-stratifications.tsv",
        fasta = FASTA
    output:
        out = directory(outdir + "/stratification/{caller}/{name}")
    log:
        outdir + "/stratification/{caller}/{name}.log"
    threads:
        8
    shell:
        """(
        mkdir -p {output.out}
        set +u; source activate happy
        hap.py -r {input.fasta} -o {output.out}/benchmark --threads {threads} \
            --engine xcmp --no-roc --stratification {input.tsv} \
            {input.vcf1} {input.vcf2} ) &> {log}
        """
