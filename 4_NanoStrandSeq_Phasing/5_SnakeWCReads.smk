#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = assembly_dir + "/wc"


rule all:
    input:
        expand(outdir + "/{cell}", cell=cells),

rule fetch_wc_regions:
    input:
        bam = assembly_dir + "/prepare/bams/{cell}.bam",
        bed = assembly_dir + "/blackhole/blacklist.bed.gz"
    output:
        out = directory(outdir + "/{cell}")
    log:
        outdir + "/{cell}.log"
    shell:
        """
        sstools FetchWCRegion -b {input.bed} {input.bam} {output.out} &> {log}
        """


