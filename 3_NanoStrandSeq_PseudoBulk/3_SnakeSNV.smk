#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
COVS1 = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]
tmp = []
for name in NAMES:
    if name.endswith(".full"):
        tmp.append(name)
    else:
        if int(name.split(".")[1].split("-")[0][3:]) in COVS1:
            tmp.append(name)
NAMES = tmp
NAMES = list(filter(lambda x: "RmDupFlag" not in x, NAMES))
CHROMS = ["chr%d" % c for c in range(1, 23)]
CALLERS = ["clair2", "longshot", "nanocaller"][2:]
ENGINES = ["xcmp"]

rule all:
    input:
        # expand(OUTDIR + "/snvs/chroms/{caller}/{name}/{chrom}.vcf.gz", caller=CALLERS, name=NAMES, chrom=CHROMS),
        expand(OUTDIR + "/snvs/chroms/{caller}/{name}/chrX.vcf.gz", caller=["nanocaller"], name=["PacBio.full"]), # For chrX only.
        expand(OUTDIR + "/snvs/concated/{caller}/{name}.vcf.gz", caller=CALLERS, name=NAMES),
        expand(OUTDIR + "/snvs/benchmark/{caller}/{name}.json", caller=CALLERS, name=NAMES),
        # expand(OUTDIR + "/snvs/benchmark_{engine}/{caller}/{name}", engine=ENGINES, caller=CALLERS, name=NAMES),

## Clair2

def get_clair2_model(name):
    if name.split(".")[0] == "PacBio":
        return CLAIR_MODEL_CCS
    else:
        return CLAIR_MODEL_ONT

rule clair2:
    input:
        fasta = FASTA,
        bam = OUTDIR + "/bams/{name}.bam"
    output:
        tmp = temp(OUTDIR + "/snvs/chroms/clair2/{name}/{chrom}.vcf"),
        vcf = OUTDIR + "/snvs/chroms/clair2/{name}/{chrom}.vcf.gz"
    log:
        OUTDIR + "/snvs/chroms/clair2/{name}/{chrom}.log"
    conda:
        "clair"
    params:
        model = lambda wildcards: get_clair2_model(wildcards.name)
    threads:
        8
    shell:
        """(
        clair.py callVarBam \
            --chkpnt_fn {params.model} \
            --ref_fn {input.fasta} \
            --bam_fn {input.bam} \
            --ctgName {wildcards.chrom} \
            --sampleName clair2 \
            --threads {threads} \
            --call_fn {output.tmp}
        bgzip -c {output.tmp} > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

## NanoCaller

def get_nanocaller_parameters(name):
    name = name.split(".")[0]
    if name == "PacBio":
        return "-p ccs --snp_model CCS-HG001 --indel_model CCS-HG001"
    if name == "Ultralong":
        return "-p ul_ont --snp_model ONT-HG001 --indel_model ONT-HG001"
    return "-p ont --snp_model ONT-HG001 --indel_model ONT-HG001"

rule nanocaller: # do not support depth <= 1
    input:
        fasta = FASTA,
        bam = OUTDIR + "/bams/{name}.bam"
    output:
        out = temp(directory(OUTDIR + "/snvs/chroms/nanocaller/{name}/{chrom}")),
        vcf = OUTDIR + "/snvs/chroms/nanocaller/{name}/{chrom}.vcf.gz"
    log:
        OUTDIR + "/snvs/chroms/nanocaller/{name}/{chrom}.log"
    conda:
        "NanoCaller"
    params:
        sup = lambda wildcards: get_nanocaller_parameters(wildcards.name)
    threads:
        12
    shell:
        """(
        python ${{HOME}}/software/NanoCaller-master/scripts/NanoCaller.py \
            -bam {input.bam} -o {output.out} -chrom {wildcards.chrom} \
            -ref {input.fasta} -cpu {threads} {params.sup} -sample nanocaller
        cp {output.out}/variant_calls.final.vcf.gz {output.vcf}
        cp {output.out}/variant_calls.final.vcf.gz.tbi {output.vcf}.tbi ) &> {log}
        """

## Longshot

rule longshot:
    input:
        fasta = FASTA,
        bam = OUTDIR + "/bams/{name}.bam"
    output:
        tmp = temp(OUTDIR + "/snvs/chroms/longshot/{name}/{chrom}.vcf"),
        vcf = OUTDIR + "/snvs/chroms/longshot/{name}/{chrom}.vcf.gz"
    log:
        OUTDIR + "/snvs/chroms/longshot/{name}/{chrom}.log"
    conda:
        "longshot"
    threads:
        4
    shell:
        """(
        longshot -A -r {wildcards.chrom} --sample_id longshot --bam {input.bam} \
            --ref {input.fasta} --out {output.tmp}
        bgzip -c {output.tmp} > {output.vcf} 
        tabix -p vcf {output.vcf} ) &> {log}
        """

## Merge chrom.vcf.gz

rule concat_chrom_vcfs:
    input:
        vcfs = lambda wildcards: [OUTDIR + "/snvs/chroms/{caller}/{name}/%s.vcf.gz" % c for c in CHROMS]
    output:
        vcf = OUTDIR + "/snvs/concated/{caller}/{name}.vcf.gz"
    log:
        OUTDIR + "/snvs/concated/{caller}/{name}.log"
    shell:
        """(
        bcftools concat -a {input.vcfs} | bcftools sort | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

# Benchmark SNPs

rule benchmark_snp:
    input:
        vcf1 = SNP_VCF,
        vcf2 = OUTDIR + "/snvs/concated/{caller}/{name}.vcf.gz",
        bed = SNP_BED
    output:
        txt = OUTDIR + "/snvs/benchmark/{caller}/{name}.json"
    log:
        OUTDIR + "/snvs/benchmark/{caller}/{name}.log"
    threads:
        8
    shell:
        """
        sstools BenchmarkSNP -t {threads} -o {output.txt} -b {input.bed} {input.vcf1} {input.vcf2}
        """

# Benchmark SNPs by hap.py

rule happy:
    input:
        vcf1 = SNP_VCF,
        vcf2 = OUTDIR + "/snvs/concated/{caller}/{name}.vcf.gz",
        bed = SNP_BED,
        fasta = FASTA
    output:
        out = directory(OUTDIR + "/snvs/happy/{caller}/{name}")
    log:
        OUTDIR + "/snvs/happy/{caller}/{name}.log"
    conda:
        "happy"
    threads:
        8
    shell:
        """(
        mkdir -p {output.out}
        hap.py -r {input.fasta} -o {output.out}/benchmark \
            --threads {threads} --engine xcmp -T {input.bed} \
            {input.vcf1} {input.vcf2} ) &> {log}
        """

rule benchmark_xcmp:
    input:
        vcf1 = SNP_VCF,
        vcf2 = OUTDIR + "/snvs/concated/{caller}/{name}.vcf.gz",
        bed = SNP_BED,
        fa = FASTA
    output:
        out = directory(OUTDIR + "/snvs/benchmark_{engine}/{caller}/{name}")
    log:
        OUTDIR + "/snvs/benchmark_{engine}/{caller}/{name}.log"
    conda:
        "happy"
    threads:
        8
    shell:
        """(
        mkdir -p {output.out}
        hap.py -r {input.fa} \
            -o {output.out}/benchmark \
            --threads {threads} \
            --engine {wildcards.engine} \
            -T {input.bed} \
            {input.vcf1} {input.vcf2} ) &> {log}
        """