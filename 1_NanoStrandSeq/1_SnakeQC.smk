#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "data/datasets"
outdir = "results/qc"

rule all:
    input:
        expand(outdir + "/fastqc/{run}_fastqc.html", run=runs),
        expand(outdir + "/lengths/{run}.txt", run=runs),
        expand(outdir + "/report/{run}.tsv", run=runs),

rule fastqc:
    input:
        fq = indir + "/{run}.fastq.gz"
    output:
        fq = temp(outdir + "/fastqc/{run}.fastq"),
        html = outdir + "/fastqc/{run}_fastqc.html"
    shell:
        """
        ./scripts/qc/get_partial_reads.py {input.fq} {output.fq}
        fastqc {output.fq} > /dev/null 2>&1
        """  

rule stat_read_length:
    input:
        fq = indir + "/{run}.fastq.gz"
    output:
        txt = outdir + "/lengths/{run}.txt"
    shell:
        """
        ./scripts/qc/stat_read_length.py {input.fq} > {output.txt}
        """

rule report_length_summary:
    input:
        tsv = rules.stat_read_length.output.txt
    output:
        tsv = outdir + "/report/{run}.tsv"
    shell:
        """
        ./scripts/qc/report_length_summary.py {input.tsv} > {output.tsv}
        """
