#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"

COVS1 = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]
tmp = []
for name in NAMES:
    if name.endswith(".full"):
        tmp.append(name)
    else:
        if int(name.split(".")[-1].split("-")[0][3:]) in COVS1:
            tmp.append(name)
NAMES = tmp
NAMES = list(filter(lambda x: "RmDupFlag" not in x, NAMES))

rule all:
    input:
        expand(OUTDIR + "/sv/cutesv/{name}.vcf.gz", name=NAMES),
        expand(OUTDIR + "/sv/filtered/{name}.vcf.gz", name=NAMES),
        expand(OUTDIR + "/sv/quantify/{name}.tsv.gz", name=NAMES),
        expand(OUTDIR + "/sv/quantify_lite/{name}.tsv", name=NAMES),
        # expand(OUTDIR + "/sv/benchmark/{name1}_vs_{name2}.json", name1=["PacBio.full", "Ultralong.full"], name2=NAMES),

# cuteSV

def get_cutesv_parameters(name):
    if name.split(".")[0] == "PacBio":
        return "--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5"
    else:
        return "--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3"

rule cutesv:
    input:
        fasta = FASTA,
        bam = OUTDIR + "/bams/{name}.bam"
    output:
        wd = temp(directory(OUTDIR + "/sv/cutesv/{name}.wd")),
        tmp = temp(OUTDIR + "/sv/cutesv/{name}.vcf"),
        vcf = OUTDIR + "/sv/cutesv/{name}.vcf.gz"
    log:
        OUTDIR + "/sv/cutesv/{name}.log"
    conda:
        "cutesv"
    params:
        sup = lambda wildcards: get_cutesv_parameters(wildcards.name)
    threads:
        THREADS
    shell:
        """(
        mkdir -p {output.wd}
        cuteSV {params.sup} --threads {threads} \
            --sample `basename {output.vcf} .vcf.gz` \
            --retain_work_dir --report_readid \
            --min_support 1 --min_size 50 --genotype \
            {input.bam} {input.fasta} {output.tmp} {output.wd} 
        bgzip -c {output.tmp} > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule filter_vcf:
    input:
        vcf = rules.cutesv.output.vcf
    output:
        vcf = OUTDIR + "/sv/filtered/{name}.vcf.gz"
    shell:
        """
        zcat {input.vcf} | ./scripts/filter_cutesv_vcf.py | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# Quantify SVs

rule quantify_sv:
    input:
        vcf = rules.filter_vcf.output.vcf,
        bam = OUTDIR + "/bams/{name}.bam"
    output:
        tmp = temp(OUTDIR + "/sv/quantify/{name}.tsv"),
        txt = OUTDIR + "/sv/quantify/{name}.tsv.gz"
    threads:
        12
    shell:
        """
        ./scripts/quantify_sv.py -t {threads} -o {output.tmp} {input.vcf} {input.bam}
        pigz -p {threads} -c {output.tmp} > {output.txt}
        """

rule simplify_quant_tsv:
    input:
        txt = rules.quantify_sv.output.txt
    output:
        txt = OUTDIR + "/sv/quantify_lite/{name}.tsv"
    shell:
        """
        gzip -d -c {input.txt} \
            | awk -v OFS='\\t' '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}' > {output.txt}
        """

# benchmark of sv

# rule benchmark_sv:
#     input:
#         vcf1 = OUTDIR + "/sv/filtered/{name1}.vcf.gz",
#         vcf2 = OUTDIR + "/sv/filtered/{name2}.vcf.gz",
#         tsv1 = OUTDIR + "/sv/quantify_lite/{name1}.tsv",
#         tsv2 = OUTDIR + "/sv/quantify_lite/{name2}.tsv"
#     output:
#         txt = OUTDIR + "/sv/benchmark/{name1}_vs_{name2}.json"
#     log:
#         OUTDIR + "/sv/benchmark/{name1}_vs_{name2}.log"
#     params:
#         bed = "../0_IntergratedAnalysis/2_DownsampleBenchmark/data/benchmark_sv_blacklist.bed",
#         min_query_cell = lambda wildcards: 3 if "NSS" in wildcards.name2 else 0
#     shell:
#         """
#         ./scripts/benchmark_sv.py {input} {params.bed} {params.min_query_cell} {output} &> {log}
#         """