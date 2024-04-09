#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/final"
OUTDIR = "results/stat"
# RUN_CELLS = RUN_CELLS[:2]

rule all:
    input:
        expand(OUTDIR + "/preseq/{run_cell}", run_cell=RUN_CELLS)

rule preseq: # Need all reads, but mark duplicate is unnecessary
    input:
        bam = INDIR + "/{run}/{cell}.bam"
    output:
        out = directory(OUTDIR + "/preseq/{run}/{cell}")
    log:
        OUTDIR + "/preseq/{run}/{cell}.log"
    conda:
        "preseq"
    params:
        mr1 = OUTDIR + "/preseq/{run}/{cell}/{cell}.mr",
        mr2 = OUTDIR + "/preseq/{run}/{cell}/{cell}.filtered.mr",
        tsv = OUTDIR + "/preseq/{run}/{cell}/{cell}.tsv"
    shell:
        """(
        set +e; mkdir {output.out}
        to-mr -o {params.mr1} {input.bam}
        awk '$1~/^chr[0-9]+$/' {params.mr1} > {params.mr2}
        preseq gc_extrap -w 100000 -e 100000000000 -s 100000000 -o {params.tsv} {params.mr2}
        rm -f {params.mr1} {params.mr2}; exit 0 ) &> {log}
        """
