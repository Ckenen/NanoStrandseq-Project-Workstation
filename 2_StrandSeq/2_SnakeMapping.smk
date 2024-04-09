#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/prepare/cutadapt"
OUTDIR = "results/mapping"
# RUN_CELLS = RUN_CELLS[:2]

rule all:
    input:
        expand(OUTDIR + "/bwa/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/bwa/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/filtered/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/mark_region/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/mark_haplotype/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/mark_duplicate/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/mark_duplicate/{run_cell}.flagstat", run_cell=RUN_CELLS),
        expand(OUTDIR + "/remove_duplicate/{run_cell}.bam", run_cell=RUN_CELLS),
        expand(OUTDIR + "/remove_duplicate/{run_cell}.flagstat", run_cell=RUN_CELLS),

rule bwa_mem:
    input:
        fq1 = INDIR + "/{run}/{cell}_1.fastq.gz",
        fq2 = INDIR + "/{run}/{cell}_2.fastq.gz",
        idx = config["BWA_INDEX"]
    output:
        bam = OUTDIR + "/bwa/{run}/{cell}.bam"
    log:
        OUTDIR + "/bwa/{run}/{cell}.log"
    conda:
        "bwa"
    params:
        rg = '@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}'
    threads:
        8
    shell:
        """(
        bwa mem -t {threads} -R \'{params.rg}\' {input.idx}/ref {input.fq1} {input.fq2} \
            | samtools view -F 4 -u - \
            | samtools sort -@ {threads} -T {output.bam}_TMP - > {output.bam} 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.bwa_mem.output.bam
    output:
        bam = OUTDIR + "/filtered/{run}/{cell}.bam"
    log:
        OUTDIR + "/filtered/{run}/{cell}.log"
    threads:
        1
    shell:
        """
        sstools FilterBam -p -n '^chr([0-9]+|[XY])$' -q 20 {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_region: # optional
    input:
        bam = rules.filter_bam.output.bam,
        bed = config["SNP_BED"]
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
        vcf = config["SNP_VCF"]
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
        bam = OUTDIR + "/mark_duplicate/{run}/{cell}.bam",
        tmpdir = temp(directory(OUTDIR + "/mark_duplicate/{run}/{cell}.tmp"))
    log:
        OUTDIR + "/mark_duplicate/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        mkdir -p {output.tmpdir}
        sambamba markdup -t {threads} --tmpdir={output.tmpdir} {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule remove_duplicate:
    input:
        bam = rules.mark_duplicate.output.bam
    output:
        bam = OUTDIR + "/remove_duplicate/{run}/{cell}.bam"
    log:
        OUTDIR + "/remove_duplicate/{run}/{cell}.log"
    threads:
        4
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
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """