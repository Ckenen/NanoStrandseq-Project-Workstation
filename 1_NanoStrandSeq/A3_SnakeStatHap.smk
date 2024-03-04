#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/mark_parental"
outdir  = "results/stat_hap"


rule all:
    input:
        expand(outdir + "/{run_cell}", run_cell=run_cells),
        # expand(outdir + "/pileup/{run_cell}", run_cell=run_cells),


rule stat_cell_haplotype:
    input:
        bam = indir + "/{run}/{cell}.bam",
        bed = lambda wildcards: get_parental_het_snv_bed(wildcards.cell),
        fasta = lambda wildcards: get_genome_fasta(wildcards.cell)
    output:
        out = directory(outdir + "/{run}/{cell}")
    log:
        outdir + "/{run}/{cell}.log"
    threads:
        8
    shell:
        """
        sstools StatCellHap -p {threads} -g {input.fasta} -a {input.bed} {input.bam} {output.out} &> {log}
        """

rule pileup:
    input:
        bam = indir + "/{run}/{cell}.bam",
        fasta = lambda wildcards: get_genome_fasta(wildcards.cell)
    output:
        out = directory(outdir + "/pileup/{run}/{cell}")
    log:
        outdir + "/pileup/{run}/{cell}.log"
    threads:
        threads
    shell:
        """
        sstools Pileup -f {input.fasta} -p {threads} {input.bam} {output.out} &> {log}
        """