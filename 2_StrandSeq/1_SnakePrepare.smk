#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
RS = ["1", "2"]
OUTDIR = "results/prepare"

rule all:
    input:
        expand(OUTDIR + "/cutadapt/{run_cell}_{r}.fastq.gz", run_cell=RUN_CELLS, r=RS),

rule prefetch:
    output:
        sra = OUTDIR + "/sra/{run}/{srr}.sra"
    log:
        OUTDIR + "/sra/{run}/{srr}.log"
    shell:
        """
        prefetch -o {output.sra} {wildcards.srr}
        """

rule sra2fastq:
    input:
        sra = rules.prefetch.output.sra
    output:
        tmp1 = temp(OUTDIR + "/fastq/{run}/{srr}_1.fastq"),
        tmp2 = temp(OUTDIR + "/fastq/{run}/{srr}_2.fastq"),
        fq1 = OUTDIR + "/fastq/{run}/{srr}_1.fastq.gz",
        fq2 = OUTDIR + "/fastq/{run}/{srr}_2.fastq.gz"
    log:
        OUTDIR + "/fastq/{run}/{srr}.log"
    threads:
        6
    shell:
        """(
        fasterq-dump --threads {threads} --split-3 --OUTDIR `dirname {output.fq1}` {input.sra}
        pigz -p {threads} -c {output.tmp1} > {output.fq1}
        pigz -p {threads} -c {output.tmp2} > {output.fq2} ) &> {log}
        """

def get_fastqs(wildcards):
    run, cell = wildcards.run, wildcards.cell
    srr = DAT[DAT["Cell"] == cell]["SRR"].values[0]
    fq1 = OUTDIR + "/fastq/%s/%s_1.fastq.gz" % (run, srr)
    fq2 = OUTDIR + "/fastq/%s/%s_2.fastq.gz" % (run, srr)
    return fq1, fq2

rule cutadapt:
    input:
        fqs = lambda wildcards: get_fastqs(wildcards)
    output:
        fq1 = temp(OUTDIR + "/cutadapt/{run}/{cell}_1.temp.fastq"),
        fq2 = temp(OUTDIR + "/cutadapt/{run}/{cell}_2.temp.fastq"),
        fq3 = OUTDIR + "/cutadapt/{run}/{cell}_1.fastq.gz",
        fq4 = OUTDIR + "/cutadapt/{run}/{cell}_2.fastq.gz",
        txt1 = OUTDIR + "/cutadapt/{run}/{cell}.round1.log",
        txt2 = OUTDIR + "/cutadapt/{run}/{cell}.round2.log"
    conda:
        "cutadapt"
    threads:
        8
    shell:
        """
        cutadapt -j {threads} -m 30 -q 30 \
            -a AGATCGGAAGAGC -a GATCGGAAGAGC \
            -A AGATCGGAAGAGC -A GATCGGAAGAGC \
            -o {output.fq1} -p {output.fq2} \
            {input.fqs[0]} {input.fqs[1]} &> {output.txt1}
        cutadapt -j {threads} -m 30 -q 30 \
            -a AGATCGGAAGAGC -a GATCGGAAGAGC \
            -A AGATCGGAAGAGC -A GATCGGAAGAGC \
            -o {output.fq3} -p {output.fq4} \
            {output.fq1} {output.fq2} &> {output.txt2}
        """
