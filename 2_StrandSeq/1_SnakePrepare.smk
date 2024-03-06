#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
rs = ["1", "2"]
outdir = "results/prepare"

rule all:
    input:
        expand(outdir + "/fastq_renamed/{run_cell}_{r}.fastq.gz", run_cell=run_cells, r=rs),
        expand(outdir + "/cutadapt/{run_cell}_{r}.fastq.gz", run_cell=run_cells, r=rs),

rule prefetch:
    output:
        sra = outdir + "/sra/{run}/{srr}.sra"
    log:
        outdir + "/sra/{run}/{srr}.log"
    shell:
        """
        prefetch -o {output.sra} {wildcards.srr}
        """

rule sra2fastq:
    input:
        sra = rules.prefetch.output.sra
    output:
        tmp1 = temp(outdir + "/fastq/{run}/{srr}_1.fastq"),
        tmp2 = temp(outdir + "/fastq/{run}/{srr}_2.fastq"),
        fq1 = outdir + "/fastq/{run}/{srr}_1.fastq.gz",
        fq2 = outdir + "/fastq/{run}/{srr}_2.fastq.gz"
    log:
        outdir + "/fastq/{run}/{srr}.log"
    threads:
        6
    shell:
        """(
        fasterq-dump --threads {threads} --split-3 --outdir `dirname {output.fq1}` {input.sra}
        pigz -p {threads} -c {output.tmp1} > {output.fq1}
        pigz -p {threads} -c {output.tmp2} > {output.fq2} ) &> {log}
        """

def get_fastqs(wildcards):
    run, cell = wildcards.run, wildcards.cell
    srr = dat[dat["Cell"] == cell]["SRR"].values[0]
    fq1 = outdir + "/fastq/%s/%s_1.fastq.gz" % (run, srr)
    fq2 = outdir + "/fastq/%s/%s_2.fastq.gz" % (run, srr)
    return fq1, fq2

rule link_fastqs:
    input:
        fqs = lambda wildcards: get_fastqs(wildcards)
    output:
        fq1 = outdir + "/fastq_renamed/{run}/{cell}_1.fastq.gz",
        fq2 = outdir + "/fastq_renamed/{run}/{cell}_2.fastq.gz"
    shell:
        """
        ln -s `readlink -f {input.fqs[0]}` {output.fq1}
        ln -s `readlink -f {input.fqs[1]}` {output.fq2}
        """

rule cutadapt:
    input:
        fq1 = rules.link_fastqs.output.fq1,
        fq2 = rules.link_fastqs.output.fq2
    output:
        fq1 = temp(outdir + "/cutadapt/{run}/{cell}_1.temp.fastq"),
        fq2 = temp(outdir + "/cutadapt/{run}/{cell}_2.temp.fastq"),
        fq3 = outdir + "/cutadapt/{run}/{cell}_1.fastq.gz",
        fq4 = outdir + "/cutadapt/{run}/{cell}_2.fastq.gz",
        txt1 = outdir + "/cutadapt/{run}/{cell}.round1.log",
        txt2 = outdir + "/cutadapt/{run}/{cell}.round2.log"
    threads:
        8
    shell:
        """
        cutadapt -j {threads} -m 30 -q 30 \
            -a AGATCGGAAGAGC -a GATCGGAAGAGC \
            -A AGATCGGAAGAGC -A GATCGGAAGAGC \
            -o {output.fq1} -p {output.fq2} \
            {input.fq1} {input.fq2} &> {output.txt1}
        cutadapt -j {threads} -m 30 -q 30 \
            -a AGATCGGAAGAGC -a GATCGGAAGAGC \
            -A AGATCGGAAGAGC -A GATCGGAAGAGC \
            -o {output.fq3} -p {output.fq4} \
            {output.fq1} {output.fq2} &> {output.txt2}
        """
