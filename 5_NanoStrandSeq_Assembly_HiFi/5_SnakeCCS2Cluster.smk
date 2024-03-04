#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"

clusters = [line.strip().split("\t")[0] for line in open(cluster_fasta + ".fai")]
outdir = "results/ccs2cluster"
# clusters = ["cluster_26"]

rule all:
    input:
        outdir + "/minimap2/cluster.map-hifi.mmi",
        expand(outdir + "/minimap2/mapped/{sample}.bam", sample=hifi_samples),
        outdir + "/minimap2/mapped_merged.bam",
        outdir + "/minimap2/mapped_merged_filtered.bam",
        expand(outdir + "/nanocaller/{cluster}.vcf.gz", cluster=clusters),
        outdir + "/nanocaller.vcf.gz",

rule build_index:
    input:
        fa = cluster_fasta
    output:
        mmi = outdir + "/minimap2/cluster.map-hifi.mmi"
    log:
        outdir + "/minimap2/cluster.map-hifi.log"
    threads:
        24
    shell:
        """
        minimap2 -t {threads} -x map-hifi -d {output.mmi} {input.fa} &> {log}
        """

rule minimap2:
    input:
        mmi = rules.build_index.output.mmi,
        fq = lambda wildcards: get_ccs_fastq(wildcards.sample)
    output:
        bam = outdir + "/minimap2/mapped/{sample}.bam"
    log:
        outdir + "/minimap2/mapped/{sample}.log"
    params:
        rg = "@RG\\tID:{sample}\\tLB:{sample}\\tSM:{sample}"
    threads:
        12
    shell:
        """(
        minimap2 -ax map-ont --MD -t {threads} -R '{params.rg}' {input.mmi} {input.fq} \
            | samtools view -@ {threads} -u -F 4 - \
            | samtools sort -@ {threads} -T {output.bam}_TMP -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule merge_bams:
    input:
        bams = expand(outdir + "/minimap2/mapped/{sample}.bam", sample=hifi_samples)
    output:
        bam = outdir + "/minimap2/mapped_merged.bam"
    threads:
        8
    shell:
        """
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        """

rule filter_bam:
    input:
        bam = rules.merge_bams.output.bam
    output:
        bam = outdir + "/minimap2/mapped_merged_filtered.bam"
    log:
        outdir + "/minimap2/mapped_merged_filtered.log"
    threads:
        8
    shell:
        """
        samtools view -@ {threads} -F 2308 -q 30 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule nanocaller:
    input:
        fasta = cluster_fasta,
        bam = rules.filter_bam.output.bam
    output:
        out = temp(directory(outdir + "/nanocaller/{cluster}")),
        vcf = outdir + "/nanocaller/{cluster}.vcf.gz"
    log:
        outdir + "/nanocaller/{cluster}.log"
    threads:
        12
    shell:
        """(
        set +u; source activate NanoCaller
        python /home/chenzonggui/software/NanoCaller-master/scripts/NanoCaller.py \
            -bam {input.bam} -o {output.out} \
            -chrom {wildcards.cluster} \
            -ref {input.fasta} -cpu {threads} \
            -p ccs --snp_model CCS-HG001 --indel_model CCS-HG001 \
            -sample nanocaller
        conda deactivate
        cp {output.out}/variant_calls.final.vcf.gz {output.vcf}
        cp {output.out}/variant_calls.final.vcf.gz.tbi {output.vcf}.tbi ) &> {log}
        """

rule merge_vcf:
    input:
        vcfs = expand(outdir + "/nanocaller/{cluster}.vcf.gz", cluster=clusters)
    output:
        vcf = outdir + "/nanocaller.vcf.gz"
    shell:
        """
        bcftools concat -a {input.vcfs} | bcftools sort | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# rule whatshap_phase:
#     input:
#         fa = cluster_fasta,
#         vcf = rules.nanocaller.output.vcf,
#         bam = rules.filter_bam.output.bam
#     output:
#         vcf = outdir + "/nanocaller_local_phased/{cluster}.vcf.gz"
#     log:
#         outdir + "/nanocaller_local_phased/{cluster}.log"
#     shell:
#         """
#         whatshap phase --ignore-read-groups --chromosome {wildcards.cluster} -r {input.fa} {input.vcf} {input.bam} | bgzip -  &> {log}
#         """
