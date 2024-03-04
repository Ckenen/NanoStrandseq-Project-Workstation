#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/minimap2"
outdir = "results/ashleys"
# run_cells = run_cells[:10]

rule all:
    input:
        expand(outdir + "/features/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/predict/{run_cell}_default.tsv", run_cell=run_cells),

rule features:
    input:
        bam = "results/mapping/minimap2/{run}/{cell}.bam"
    output:
        tsv = outdir + "/features/{run}/{cell}.tsv"
    log:
        outdir + "/features/{run}/{cell}.log"
    shell:
        """
        set +u; source activate ashleys
        ashleys.py features -f {input.bam} -w 5000000 2000000 1000000 800000 600000 400000 200000 -o {output.tsv} &> {log}
        conda deactivate
        """

rule predict:
    input:
        tsv = rules.features.output.tsv
    output:
        tsv1 = outdir + "/predict/{run}/{cell}_default.tsv",
        tsv2 = outdir + "/predict/{run}/{cell}_stringent.tsv"
    log:
        outdir + "/predict/{run}/{cell}.log"
    params:
        model_default = "/home/chenzonggui/software/ashleys-qc-pipeline-main/workflow/ashleys_models/svc_default.pkl",
        model_stringent = "/home/chenzonggui/software/ashleys-qc-pipeline-main/workflow/ashleys_models/svc_stringent.pkl"
    shell:
        """(
        set +u; source activate ashleys
        ashleys.py predict -p {input.tsv} -o {output.tsv1} -m {params.model_default}
        ashleys.py predict -p {input.tsv} -o {output.tsv2} -m {params.model_stringent}
        conda deactivate ) &> {log}
        """
