#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/prepare/cutadapt"
outdir = "results/mapping"
# run_cells = run_cells[:1]

rule all:
    input:
        expand(outdir + "/bwa/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/bwa/{run_cell}.flagstat", run_cell=run_cells),
        # expand(outdir + "/filtered/{run_cell}.bam", run_cell=run_cells),
        # expand(outdir + "/mark_region/{run_cell}.bam", run_cell=run_cells),
        # expand(outdir + "/mark_haplotype/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/mark_duplicate/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/mark_duplicate/{run_cell}.flagstat", run_cell=run_cells),
        expand(outdir + "/remove_duplicate/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/remove_duplicate/{run_cell}.flagstat", run_cell=run_cells),

rule bwa_mem:
    input:
        fq1 = indir + "/{run}/{cell}_1.fastq.gz",
        fq2 = indir + "/{run}/{cell}_2.fastq.gz",
        idx = config["bwa_index"]
    output:
        bam = outdir + "/bwa/{run}/{cell}.bam"
    log:
        outdir + "/bwa/{run}/{cell}.log"
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
        bam = outdir + "/filtered/{run}/{cell}.bam"
    log:
        outdir + "/filtered/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools FilterBam -p -n '^chr([0-9]+|[XY])$' -q 20 {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_region: # optional
    input:
        bam = rules.filter_bam.output.bam,
        bed = config["benchmark_bed"]
    output:
        bam = outdir + "/mark_region/{run}/{cell}.bam"
    log:
        outdir + "/mark_region/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools MarkRegion -n XH {input.bam} {input.bed} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_haplotype: # optional
    input:
        bam = rules.mark_region.output.bam,
        vcf = config["benchmark_vcf"]
    output:
        bam = outdir + "/mark_haplotype/{run}/{cell}.bam"
    log:
        outdir + "/mark_haplotype/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools MarkHaplotype -n XP {input.bam} {input.vcf} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_duplicate:
    input:
        bam = rules.mark_haplotype.output.bam
    output:
        bam = outdir + "/mark_duplicate/{run}/{cell}.bam",
        tmpdir = temp(directory(outdir + "/mark_duplicate/{run}/{cell}.tmp"))
    log:
        outdir + "/mark_duplicate/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        mkdir -p {output.tmpdir}
        sambamba markdup -t {threads} --tmpdir={output.tmpdir} {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

# DEPRECATED (too slow)

# rule mark_duplicate:
#     input:
#         bam = rules.filter_bam.output.bam
#     output:
#         bam = outdir + "/mark_duplicate/{run}/{cell}.bam",
#         txt = outdir + "/mark_duplicate/{run}/{cell}.metrics"
#     log:
#         outdir + "/mark_duplicate/{run}/{cell}.log"
#     threads:
#         4
#     shell:
#         """
#         picard MarkDuplicates --REMOVE_DUPLICATES false -I {input.bam} -M {output.txt} -O {output.bam} &> {log}
#         samtools index -@ {threads} {output.bam}
#         """

rule remove_duplicate:
    input:
        bam = rules.mark_duplicate.output.bam
    output:
        bam = outdir + "/remove_duplicate/{run}/{cell}.bam"
    log:
        outdir + "/remove_duplicate/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} -F 1024 -b {input.bam} > {output.bam}
        samtools index -@ {threads} {output.bam}
        """

# Common rule

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