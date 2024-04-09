#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "data/datasets"
OUTDIR = "results/qc"

rule all:
    input:
        expand(OUTDIR + "/fastqc/{run}_fastqc.html", run=RUNS),
        expand(OUTDIR + "/lengths/{run}.txt", run=RUNS),
        expand(OUTDIR + "/report/{run}.tsv", run=RUNS),

rule fastqc:
    input:
        fq = INDIR + "/{run}.fastq.gz"
    output:
        fq = temp(OUTDIR + "/fastqc/{run}.fastq"),
        html = OUTDIR + "/fastqc/{run}_fastqc.html"
    shell:
        """
        ./scripts/qc/get_partial_reads.py {input.fq} {output.fq}
        fastqc {output.fq} > /dev/null 2>&1
        """  

rule stat_read_length:
    input:
        fq = INDIR + "/{run}.fastq.gz"
    output:
        txt = OUTDIR + "/lengths/{run}.txt"
    shell:
        """
        ./scripts/qc/stat_read_length.py {input.fq} > {output.txt}
        """

rule report_length_summary:
    input:
        tsv = rules.stat_read_length.output.txt
    output:
        tsv = OUTDIR + "/report/{run}.tsv"
    shell:
        """
        ./scripts/qc/report_length_summary.py {input.tsv} > {output.tsv}
        """
