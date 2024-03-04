#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/mark_parental"
outdir = "results/stat"

rule all:
    input:
        expand(outdir + "/preseq/{run_cell}", run_cell=run_cells)


rule preseq: # Need all reads, but mark duplicate is unnecessary
    input:
        bam = indir + "/{run}/{cell}.bam"
    output:
        out = directory(outdir + "/preseq/{run}/{cell}")
    log:
        outdir + "/preseq/{run}/{cell}.log"
    params:
        mr1 = outdir + "/preseq/{run}/{cell}/{cell}.mr",
        mr2 = outdir + "/preseq/{run}/{cell}/{cell}.filtered.mr",
        tsv = outdir + "/preseq/{run}/{cell}/{cell}.tsv"
    shell:
        """(
        set +e; mkdir {output.out}
        to-mr -o {params.mr1} {input.bam}
        awk '$1~/^chr[0-9]+$/' {params.mr1} > {params.mr2}
        preseq gc_extrap -w 100000 -e 100000000000 -s 100000000 -o {params.tsv} {params.mr2}
        rm {params.mr1} {params.mr2}; exit 0 ) &> {log}
        """
