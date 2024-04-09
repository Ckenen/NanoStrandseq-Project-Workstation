#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/final"
OUTDIR = "results/counts"
# RUN_CELLS = RUN_CELLS[:1]

rule all:
    input:
        expand(OUTDIR + "/stat_bin_reads/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/plot_bin_reads/{run_cell}.pdf", run_cell=RUN_CELLS)

rule stat_bin_reads:
    input:
        bam = INDIR + "/{run}/{cell}.bam"
    output:
        tsv = OUTDIR + "/stat_bin_reads/{run}/{cell}.tsv"
    log:
        OUTDIR + "/stat_bin_reads/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        sstools StatBinRead --width 1000000 --remove-duplicates \
            --threads {threads} {input.bam} {output.tsv} &> {log}
        """

rule plot_bin_reads:
    input:
        tsv = rules.stat_bin_reads.output.tsv
    output:
        pdf1 = OUTDIR + "/plot_bin_reads/{run}/{cell}.pdf",
        pdf2 = OUTDIR + "/plot_bin_reads/{run}/{cell}.trimmed.pdf"
    log:
        OUTDIR + "/plot_bin_reads/{run}/{cell}.log"
    shell:
        """(
        sstools PlotBinRead {input.tsv} {output.pdf1}
        sstools PlotBinRead --trim {input.tsv} {output.pdf2} ) &> {log}
        """
