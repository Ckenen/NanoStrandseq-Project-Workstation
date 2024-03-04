#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"

covs1 = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]
tmp = []
for name in names:
    if name.endswith(".full"):
        tmp.append(name)
    else:
        if int(name.split(".")[1].split("-")[0][3:]) in covs1:
            tmp.append(name)

names = tmp

chroms = ["chr%d" % c for c in range(1, 23)]
callers = ["clair2", "longshot", "nanocaller"][2:]
engines = ["xcmp"]

rule all:
    input:
        # expand(outdir + "/snvs/chroms/{caller}/{name}/{chrom}.vcf.gz", caller=callers, name=names, chrom=chroms),
        expand(outdir + "/snvs/chroms/{caller}/{name}/chrX.vcf.gz", caller=["nanocaller"], name=["PacBio.full"]), # For chrX only.
        expand(outdir + "/snvs/concated/{caller}/{name}.vcf.gz", caller=callers, name=names),
        expand(outdir + "/snvs/benchmark/{caller}/{name}.json", caller=callers, name=names),
        #expand(outdir + "/snvs/benchmark_{engine}/{caller}/{name}", engine=engines, caller=callers, name=names),

## Clair2

def get_clair2_model(name):
    if name.split(".")[0] == "PacBio":
        return "/home/chenzonggui/software/princess/bin/modules/ccs/model"
    else:
        return "/home/chenzonggui/software/princess/bin/modules/ont/model"

rule clair2:
    input:
        fasta = FASTA,
        bam = outdir + "/bams/{name}.bam"
    output:
        tmp = temp(outdir + "/snvs/chroms/clair2/{name}/{chrom}.vcf"),
        vcf = outdir + "/snvs/chroms/clair2/{name}/{chrom}.vcf.gz"
    log:
        outdir + "/snvs/chroms/clair2/{name}/{chrom}.log"
    params:
        model = lambda wildcards: get_clair2_model(wildcards.name)
    threads:
        8
    shell:
        """(
        set +u; source activate clair
        clair.py callVarBam \
            --chkpnt_fn {params.model} \
            --ref_fn {input.fasta} \
            --bam_fn {input.bam} \
            --ctgName {wildcards.chrom} \
            --sampleName clair2 \
            --threads {threads} \
            --call_fn {output.tmp}
        conda deactivate
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
        bam = outdir + "/bams/{name}.bam"
    output:
        out = temp(directory(outdir + "/snvs/chroms/nanocaller/{name}/{chrom}")),
        vcf = outdir + "/snvs/chroms/nanocaller/{name}/{chrom}.vcf.gz"
    log:
        outdir + "/snvs/chroms/nanocaller/{name}/{chrom}.log"
    params:
        sup = lambda wildcards: get_nanocaller_parameters(wildcards.name)
    threads:
        12
    shell:
        """(
        set +u; source activate NanoCaller
        python /home/chenzonggui/software/NanoCaller-master/scripts/NanoCaller.py \
            -bam {input.bam} -o {output.out} -chrom {wildcards.chrom} \
            -ref {input.fasta} -cpu {threads} {params.sup} -sample nanocaller
        conda deactivate
        cp {output.out}/variant_calls.final.vcf.gz {output.vcf}
        cp {output.out}/variant_calls.final.vcf.gz.tbi {output.vcf}.tbi ) &> {log}
        """

## Longshot

rule longshot:
    input:
        fasta = FASTA,
        bam = outdir + "/bams/{name}.bam"
    output:
        tmp = temp(outdir + "/snvs/chroms/longshot/{name}/{chrom}.vcf"),
        vcf = outdir + "/snvs/chroms/longshot/{name}/{chrom}.vcf.gz"
    log:
        outdir + "/snvs/chroms/longshot/{name}/{chrom}.log"
    threads:
        4
    shell:
        """(
        set +u; source activate longshot
        longshot -A -r {wildcards.chrom} --sample_id longshot --bam {input.bam} \
            --ref {input.fasta} --out {output.tmp}
        conda deactivate 
        bgzip -c {output.tmp} > {output.vcf} 
        tabix -p vcf {output.vcf} ) &> {log}
        """

## Merge chrom.vcf.gz

rule concat_chrom_vcfs:
    input:
        vcfs = lambda wildcards: [outdir + "/snvs/chroms/{caller}/{name}/%s.vcf.gz" % c for c in chroms]
    output:
        vcf = outdir + "/snvs/concated/{caller}/{name}.vcf.gz"
    log:
        outdir + "/snvs/concated/{caller}/{name}.log"
    shell:
        """(
        bcftools concat -a {input.vcfs} | bcftools sort | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

# Benchmark SNPs

rule benchmark_snp:
    input:
        vcf1 = BENCHMARK_VCF,
        vcf2 = outdir + "/snvs/concated/{caller}/{name}.vcf.gz",
        bed = BENCHMARK_BED
    output:
        txt = outdir + "/snvs/benchmark/{caller}/{name}.json"
    log:
        outdir + "/snvs/benchmark/{caller}/{name}.log"
    threads:
        8
    shell:
        """
        sstools BenchmarkSNP -t {threads} -o {output.txt} -b {input.bed} {input.vcf1} {input.vcf2}
        """

# Benchmark SNPs by hap.py

rule happy:
    input:
        vcf1 = BENCHMARK_VCF,
        vcf2 = outdir + "/snvs/concated/{caller}/{name}.vcf.gz",
        bed = BENCHMARK_BED,
        fasta = FASTA
    output:
        out = directory(outdir + "/snvs/happy/{caller}/{name}")
    log:
        outdir + "/snvs/happy/{caller}/{name}.log"
    threads:
        8
    shell:
        """(
        mkdir -p {output.out}
        set +u; source activate happy
        hap.py -r {input.fasta} -o {output.out}/benchmark \
            --threads {threads} --engine xcmp -T {input.bed} \
            {input.vcf1} {input.vcf2} ) &> {log}
        """

# rule benchmark_xcmp:
#     input:
#         vcf1 = BENCHMARK_VCF,
#         vcf2 = outdir + "/snvs/concated/{caller}/{name}.vcf.gz",
#         bed = BENCHMARK_BED,
#         fasta = FASTA
#     output:
#         out = directory(outdir + "/snvs/benchmark_{engine}/{caller}/{name}")
#     log:
#         outdir + "/snvs/benchmark_{engine}/{caller}/{name}.log"
#     threads:
#         8
#     shell:
#         """(
#         mkdir -p {output.out}
#         set +u; source activate happy
#         hap.py -r {input.fasta} \
#             -o {output.out}/benchmark \
#             --threads {threads} \
#             --engine {wildcards.engine} \
#             -T {input.bed} \
#             {input.vcf1} {input.vcf2} ) &> {log}
#         """