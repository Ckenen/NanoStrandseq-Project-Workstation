#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/mark_parental"
outdir = "results/breakpointR"
# run_cells = run_cells[:10]

rule all:
    input:
        expand(outdir + "/breakpointr/{run_cell}", run_cell=run_cells),
        expand(outdir + "/libMetrics/{run_cell}.txt", run_cell=run_cells),

rule breakpointR:
    input:
        bam = indir + "/{run}/{cell}.bam"
    output:
        out = directory(outdir + "/breakpointr/{run}/{cell}")
    log:
        outdir + "/breakpointr/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        bam=`readlink -f {input.bam}`
        name=`basename {input.bam}`
        mkdir -p {output.out}.tmp
        ln -f -s ${{bam}} {output.out}.tmp/${{name}}
        ln -f -s ${{bam}}.bai {output.out}.tmp/${{name}}.bai
        set +u; source activate breakpointR
        ./scripts/breakpointR.SE.R {output.out}.tmp {threads} {output.out} 
        rm -r {output.out}.tmp
        conda deactivate ) &> {log}
        """

rule get_lib_metrics:
    input:
        rules.breakpointR.output.out
    output:
        txt = outdir + "/libMetrics/{run}/{cell}.txt"
    log:
        outdir + "/libMetrics/{run}/{cell}.log"
    shell:
        """
        set +u; source activate breakpointR
        ./scripts/get_lib_metrics.R {input}/data/{wildcards.cell}.bam.RData {output.txt} &> {log}
        conda deactivate
        """
