#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/demux/trimmed"
OUTDIR = "results/mapping"
# RUN_CELLS = RUN_CELLS[:2]

rule all:
    input:
        expand(OUTDIR + "/minimap2/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/minimap2/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/filtered/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/mark_region/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/mark_haplotype/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/mark_duplicate/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/mark_duplicate/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/remove_duplicate/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/remove_duplicate/{run_cell}.flagstat", run_cell=RUN_CELLS),

rule minimap2:
    input:
        fq = INDIR + "/{run}/{cell}/trimmed.fastq.gz",
        mmi = lambda wildcards: get_mmi(wildcards.cell)
    output:
        bam = OUTDIR + "/minimap2/{run}/{cell}.bam"
    log:
        OUTDIR + "/minimap2/{run}/{cell}.log"
    conda:
        "minimap2"
    params:
        rg = "@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}"
    threads:
        THREADS
    shell:
        """(
        minimap2 -ax map-ont --MD -t {threads} -R '{params.rg}' 、
            --secondary=no {input.mmi} {input.fq} \
            | samtools view -@ {threads} -u -F 4 - \
            | samtools sort -@ {threads} -T {output.bam}_TMP -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = OUTDIR + "/filtered/{run}/{cell}.bam"
    log:
        OUTDIR + "/filtered/{run}/{cell}.log"
    threads:
        1
    shell:
        """(
        sstools FilterBam -n '^chr([0-9]+|[XY])$' -q 30 {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule mark_region: # optional
    input:
        bam = rules.filter_bam.output.bam,
        bed = lambda wildcards: get_snp_bed(wildcards.cell)
    output:
        bam = OUTDIR + "/mark_region/{run}/{cell}.bam"
    log:
        OUTDIR + "/mark_region/{run}/{cell}.log"
    threads:
        1
    shell:
        """
        sstools MarkRegion -n XH {input.bam} {input.bed} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_haplotype: # optional
    input:
        bam = rules.mark_region.output.bam,
        vcf = lambda wildcards: get_snp_vcf(wildcards.cell)
    output:
        bam = OUTDIR + "/mark_haplotype/{run}/{cell}.bam"
    log:
        OUTDIR + "/mark_haplotype/{run}/{cell}.log"
    threads:
        1
    shell:
        """
        sstools MarkHaplotype -n XP {input.bam} {input.vcf} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_duplicate:
    input:
        bam = rules.mark_haplotype.output.bam
    output:
        bam = OUTDIR + "/mark_duplicate/{run}/{cell}.bam"
    log:
        OUTDIR + "/mark_duplicate/{run}/{cell}.log"
    threads:
        1
    shell:
        """
        sstools MarkDuplicate -d 20 -s {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule remove_duplicate:
    input:
        bam = rules.mark_duplicate.output.bam
    output:
        bam = OUTDIR + "/remove_duplicate/{run}/{cell}.bam"
    threads:
        1
    shell:
        """
        samtools view -@ {threads} -F 1024 -b {input.bam} > {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    threads:
        1
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """
