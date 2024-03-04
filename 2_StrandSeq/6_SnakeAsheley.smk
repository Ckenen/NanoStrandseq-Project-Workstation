#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/mark_duplicates"
outdir = "results/ashleys"
models = ["default", "stringent"]

rule all:
    input:
        expand(outdir + "/features/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/predict/{run_cell}_{model}.tsv", run_cell=run_cells, model=models),

rule features:
    input:
        bam = indir + "/{run}/{cell}.bam"
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
        tsv = rules.features.output.tsv,
        model1 = "/home/chenzonggui/software/ashleys-qc-pipeline-main/workflow/ashleys_models/svc_default.pkl",
        model2 = "/home/chenzonggui/software/ashleys-qc-pipeline-main/workflow/ashleys_models/svc_stringent.pkl"
    output:
        tsv1 = outdir + "/predict/{run}/{cell}_default.tsv",
        tsv2 = outdir + "/predict/{run}/{cell}_stringent.tsv"
    log:
        outdir + "/predict/{run}/{cell}.log"
    shell:
        """(
        set +u; source activate ashleys
        ashleys.py predict -p {input.tsv} -o {output.tsv1} -m {input.model1}
        ashleys.py predict -p {input.tsv} -o {output.tsv2} -m {input.model2}
        conda deactivate ) &> {log}
        """