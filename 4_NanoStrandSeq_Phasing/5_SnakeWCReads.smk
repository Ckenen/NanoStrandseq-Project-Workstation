#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
OUTDIR = ROOT_DIR + "/wc"

rule all:
    input:
        expand(OUTDIR + "/{cell}", cell=CELLS),

rule fetch_wc_regions:
    input:
        bam = ROOT_DIR + "/prepare/bams/{cell}.bam",
        bed = ROOT_DIR + "/blackhole/blacklist.bed.gz"
    output:
        directory(OUTDIR + "/{cell}")
    log:
        OUTDIR + "/{cell}.log"
    shell:
        """
        sstools FetchWCRegion -b {input.bed} {input.bam} {output} &> {log}
        """


