#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = "results/ss2cluster"


rule all:
    input:
        # NGS
        outdir + "/ngs/bwa/cluster.index",
        expand(outdir + "/ngs/bwa/mapped/{cell}.bam", cell=strand_seq_cells),
        expand(outdir + "/ngs/filtered/{cell}.bam", cell=strand_seq_cells),
        expand(outdir + "/ngs/mark_duplicate/{cell}.bam", cell=strand_seq_cells),
        outdir + "/ngs/all_cells.bam",
        outdir + "/ngs/all_cells.bw",
        expand(outdir + "/ngs/bin_reads/{cell}.tsv", cell=strand_seq_cells),
        expand(outdir + "/ngs/bin_reads/{cell}.pdf", cell=strand_seq_cells),
        # TGS
        outdir + "/tgs/minimap2/cluster.map-ont.mmi",
        expand(outdir + "/tgs/minimap2/mapped/{cell}.bam", cell=nano_strand_seq_cells),
        expand(outdir + "/tgs/filtered/{cell}.bam", cell=nano_strand_seq_cells),
        expand(outdir + "/tgs/mark_duplicate/{cell}.bam", cell=nano_strand_seq_cells),
        outdir + "/tgs/all_cells.bam",
        outdir + "/tgs/all_cells.bw",
        expand(outdir + "/tgs/bin_reads/{cell}.tsv", cell=nano_strand_seq_cells),
        expand(outdir + "/tgs/bin_reads/{cell}.pdf", cell=nano_strand_seq_cells),
        

# NGS

rule bwa_index:
    input:
        fa = cluster_fasta
    output:
        out = directory(outdir + "/ngs/bwa/cluster.index")
    log:
        outdir + "/bwa/cluster.index.log"
    shell:
        """
        mkdir -p {output.out}
        bwa index -p {output.out}/ref {input.fa} &> {log}
        """

rule bwa_mem:
    input:
        fqs = lambda wildcards: get_strand_seq_fastqs(wildcards.cell),
        idx = rules.bwa_index.output.out
    output:
        bam = outdir + "/ngs/bwa/mapped/{cell}.bam"
    log:
        outdir + "/ngs/bwa/mapped/{cell}.log"
    params:
        rg = '@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}'
    threads:
        8
    shell:
        """(
        bwa mem -t {threads} -R \'{params.rg}\' {input.idx}/ref {input.fqs[0]} {input.fqs[1]} \
            | samtools view -F 4 -u - \
            | samtools sort -@ {threads} -T {output.bam}_TMP - > {output.bam} 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam_ngs:
    input:
        bam = rules.bwa_mem.output.bam
    output:
        bam = outdir + "/ngs/filtered/{cell}.bam"
    log:
        outdir + "/ngs/filtered/{cell}.log"
    threads:
        4
    shell:
        """(
        samtools view -@ {threads} -F 2308 -q 30 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule mark_duplicate_ngs:
    input:
        bam = rules.filter_bam_ngs.output.bam
    output:
        bam = outdir + "/ngs/mark_duplicate/{cell}.bam"
    log:
        outdir + "/ngs/mark_duplicate/{cell}.log"
    threads:
        4
    shell:
        """
        sambamba markdup -t {threads} -r {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule merge_bams_ngs:
    input:
        bams = expand(outdir + "/ngs/mark_duplicate/{cell}.bam", cell=strand_seq_cells)
    output:
        bam = outdir + "/ngs/all_cells.bam"
    threads:
        8
    shell:
        """
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam}
        """

# TGS


rule build_index:
    input:
        fa = cluster_fasta
    output:
        mmi = outdir + "/tgs/minimap2/cluster.map-ont.mmi"
    log:
        outdir + "/tgs/minimap2/cluster.map-ont.log"
    threads:
        24
    shell:
        """
        minimap2 -t {threads} -x map-ont -d {output.mmi} {input.fa} &> {log}
        """

rule minimap2:
    input:
        mmi = rules.build_index.output.mmi,
        fq = lambda wildcards: get_nano_strand_seq_fastq(wildcards.cell)
    output:
        bam = outdir + "/tgs/minimap2/mapped/{cell}.bam"
    log:
        outdir + "/tgs/minimap2/mapped/{cell}.log"
    params:
        rg = "@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}"
    threads:
        12
    shell:
        """(
        minimap2 -ax map-ont --MD -t {threads} -R '{params.rg}' {input.mmi} {input.fq} \
            | samtools view -@ {threads} -u -F 4 - \
            | samtools sort -@ {threads} -T {output.bam}_TMP -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam_tgs:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = outdir + "/tgs/filtered/{cell}.bam"
    log:
        outdir + "/tgs/filtered/{cell}.log"
    threads:
        4
    shell:
        """
        sstools FilterBam -q 60 {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_duplicate_tgs:
    input:
        bam = rules.filter_bam_tgs.output.bam,
    output:
        bam = outdir + "/tgs/mark_duplicate/{cell}.bam",
        tsv1 = outdir + "/tgs/mark_duplicate/{cell}.reads.tsv",
        tsv2 = outdir + "/tgs/mark_duplicate/{cell}.groups.tsv"
    log:
        outdir + "/tgs/mark_duplicate/{cell}.log"
    threads:
        4
    shell:
        """
        sstools MarkDuplicate --max-distance 20 --strand-sense \
            --read-matrix {output.tsv1} --group-matrix {output.tsv2} \
            {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule merge_bams_tgs:
    input:
        bams = expand(outdir + "/tgs/mark_duplicate/{cell}.bam", cell=nano_strand_seq_cells)
    output:
        bam = outdir + "/tgs/all_cells.bam"
    threads:
        8
    shell:
        """
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam}
        """

# Common

rule get_bin_reads:
    input:
        bam = outdir + "/{tech}/mark_duplicate/{cell}.bam"
    output:
        txt = outdir + "/{tech}/bin_reads/{cell}.tsv"
    shell:
        """
        ./scripts/get_bin_reads.py {input.bam} {output.txt}
        """

rule plot_bin_reads:
    input:
        txt = rules.get_bin_reads.output.txt
    output:
        pdf = outdir + "/{tech}/bin_reads/{cell}.pdf"
    shell:
        """
        ./scripts/plot_bin_reads.py {input.txt} {output.pdf}
        """

rule bam_to_bw:
    input:
        bam = "{prefix}.bam",
        txt = cluster_sizes
    output:
        td = temp(directory("{prefix}.BAM_TO_BW_SORT_TMP")),
        bed = temp("{prefix}.bed"),
        bg = temp("{prefix}.bedGraph"),
        bw = "{prefix}.bw"
    threads:
        8
    shell:
        """
        mkdir -p {output.td}
        bedtools bamtobed -i {input.bam} \
            | awk -v FS='\\t' -v OFS='\\t' '{{print $1,$2,$3}}' \
            | sort --parallel {threads} -T {output.td} -k1,1 -k2,2n -k3,3n > {output.bed}
        bedtools genomecov -bg -i {output.bed} -g {input.txt} | sort --parallel {threads} -T {output.td} -k1,1 -k2,2n > {output.bg}
        bedGraphToBigWig {output.bg} {input.txt} {output.bw}
        """