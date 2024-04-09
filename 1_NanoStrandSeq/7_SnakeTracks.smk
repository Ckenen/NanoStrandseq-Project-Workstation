#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
OUTDIR = "results/tracks"
CATEGORIES = ["all", "rmdup"][:1]
RUN_CELLS = [
    "20220708_GM12878/20220708_GM12878.sc002",
    "20220708_GM12878/20220708_GM12878.sc006",
    "20220708_GM12878/20220708_GM12878.sc010",
    "20220708_GM12878/20220708_GM12878.sc017",
    "20220708_GM12878/20220708_GM12878.sc020",
    "20220708_GM12878/20220708_GM12878.sc022",
    "20220708_GM12878/20220708_GM12878.sc029",
    "20220708_GM12878/20220708_GM12878.sc039",
    "20220708_GM12878/20220708_GM12878.sc042",
    "20220708_GM12878/20220708_GM12878.sc047",
    "20220708_GM12878/20220708_GM12878.sc055",
    "20220708_GM12878/20220708_GM12878.sc066",
    "20220708_GM12878/20220708_GM12878.sc076",
    "20220708_GM12878/20220708_GM12878.sc079",
    "20220708_GM12878/20220708_GM12878.sc087",
]

rule all:
    input:
        expand(OUTDIR + "/{cate}/{run_cell}", cate=CATEGORIES, run_cell=RUN_CELLS),

def get_bam(wildcards):
    cate = wildcards.cate
    run = wildcards.run
    cell = wildcards.cell
    if cate == "all":
        return "results/mapping/mark_duplicate/%s/%s.bam" % (run, cell)
    elif cate == "rmdup":
        return "results/mapping/remove_duplicate/%s/%s.bam" % (run, cell)

rule make_tracks:
    input:
        bam = lambda wildcards: get_bam(wildcards)
    output:
        directory(OUTDIR + "/{cate}/{run}/{cell}")
    params:
        prefix = OUTDIR + "/{cate}/{run}/{cell}/{cell}"
    shell:
        """
        mkdir {output}
        ./scripts/make_tracks.sh {input.bam} {params.prefix}
        """

