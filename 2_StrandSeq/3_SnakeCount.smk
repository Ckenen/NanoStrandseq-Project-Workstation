#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/mark_parental"
outdir = "results/counts"

rule all:
    input:
        expand(outdir + "/stat_bin_reads/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/plot_bin_reads/{run_cell}.pdf", run_cell=run_cells)

rule stat_bin_reads:
    input:
        bam = indir + "/{run}/{cell}.bam"
    output:
        tsv = outdir + "/stat_bin_reads/{run}/{cell}.tsv",
    log:
        outdir + "/stat_bin_reads/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        sstools StatBinRead --width 1000000 --remove-duplicates --threads {threads} {input.bam} {output.tsv} &> {log}
        """

rule plot_bin_reads:
    input:
        tsv = rules.stat_bin_reads.output.tsv
    output:
        pdf1 = outdir + "/plot_bin_reads/{run}/{cell}.pdf",
        pdf2 = outdir + "/plot_bin_reads/{run}/{cell}.trimmed.pdf"
    log:
        outdir + "/plot_bin_reads/{run}/{cell}.log"
    shell:
        """(
        sstools PlotBinRead {input.tsv} {output.pdf1}
        sstools PlotBinRead --trim {input.tsv} {output.pdf2} ) &> {log}
        """
