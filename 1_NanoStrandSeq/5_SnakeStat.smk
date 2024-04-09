#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "results/mapping/final"
OUTDIR = "results/stat"
# RUN_CELLS = RUN_CELLS[:2]

rule all:
    input:
        expand(OUTDIR + "/lengths/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/gc_content/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/gc_density_pkl/{species}.pkl", species=["Human", "Mouse"]),
        expand(OUTDIR + "/gc_density/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/background/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/spikiness/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/depth/{run_cell}.tsv", run_cell=RUN_CELLS),
        expand(OUTDIR + "/coverage/{run_cell}.tsv", run_cell=RUN_CELLS),

rule stat_length:
    input:
        bam = INDIR + "/{run}/{cell}.bam"
    output:
        tsv = OUTDIR + "/lengths/{run}/{cell}.tsv",
        tsv2 = OUTDIR + "/lengths/{run}/{cell}_summary.tsv"
    log:
        OUTDIR + "/lengths/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatLength -t {threads} -l SE -s {output.tsv2} {input.bam} {output.tsv} &> {log}
        """

rule stat_gc_content:
    input:
        bam = INDIR + "/{run}/{cell}.bam",
        fa = lambda wildcards: get_fasta(wildcards.cell)
    output:
        tsv = OUTDIR + "/gc_content/{run}/{cell}.tsv",
        tsv2 = OUTDIR + "/gc_content/{run}/{cell}_summary.tsv"
    log:
        OUTDIR + "/gc_content/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatGC -t {threads} -s {output.tsv2} {input.bam} {input.fa} {output.tsv} &> {log}
        """

rule prepare_genome_gc:
    input:
        fa = lambda wildcards: config["%s_FASTA" % wildcards.species.upper()]
    output:
        pkl = OUTDIR + "/gc_density_pkl/{species}.pkl"
    log:
        OUTDIR + "/gc_density_pkl/{species}.log"
    threads:
        12
    shell:
        """
        sstools StatGCDensity --prepare -t {threads} -k {output.pkl} {input.fa} &> {log}
        """

rule stat_gc_density:
    input:
        bam = INDIR + "/{run}/{cell}.bam",
        fa = lambda wildcards: get_fasta(wildcards.cell),
        pkl = lambda wildcards: OUTDIR + "/gc_density_pkl/%s.pkl" % get_species(wildcards.cell)
    output:
        tsv = OUTDIR + "/gc_density/{run}/{cell}.tsv"
    log:
        OUTDIR + "/gc_density/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatGCDensity -t {threads} -k {input.pkl} {input.bam} {input.fa} {output.tsv} &> {log}
        """

rule stat_background:
    input:
        bam = INDIR + "/{run}/{cell}.bam"
    output:
        tsv = OUTDIR + "/background/{run}/{cell}.tsv",
        tsv2 = OUTDIR + "/background/{run}/{cell}_summary.tsv"
    log:
        OUTDIR + "/background/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatBackground -t {threads} -q 60 -s {output.tsv2} {input.bam} {output.tsv} &> {log}
        """

rule stat_spikiness:
    input:
        bam = INDIR + "/{run}/{cell}.bam"
    output:
        tsv = OUTDIR + "/spikiness/{run}/{cell}.tsv"
    log:
        OUTDIR + "/spikiness/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatSpikiness -t {threads} {input.bam} {output.tsv} &> {log}
        """

rule stat_depth:
    input:
        bam = INDIR + "/{run}/{cell}.bam"
    output:
        tsv1 = OUTDIR + "/depth/{run}/{cell}.tsv",
        tsv2 = OUTDIR + "/depth/{run}/{cell}_rmdup.tsv"
    log:
        OUTDIR + "/depth/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        sstools StatDepth -t {threads} {input.bam} {output.tsv1}
        sstools StatDepth -t {threads} -r {input.bam} {output.tsv2} ) &> {log}
        """

rule stat_coverage:
    input:
        bam = INDIR + "/{run}/{cell}.bam"
    output:
        tsv = OUTDIR + "/coverage/{run}/{cell}.tsv"
    log:
        OUTDIR + "/coverage/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatCoverage -t {threads} {input.bam} {output.tsv} &> {log}
        """

# rule stat_genomic_coverage_downsample:
#     input:
#         bam = INDIR + "/{run}/{cell}.bam"
#     output:
#         tsv = OUTDIR + "/genomic_coverage_downsample/{run}/{cell}.tsv"
#     threads:
#         8
#     shell:
#         """
#         sstools StatGenomicCov -p {threads} {input.bam} {output.tsv}
#         """


# rule stat_pcr_duplicates_downsample:
#     input:
#         bam = INDIR + "/{run}/{cell}.bam",
#     output:
#         tsv = OUTDIR + "/pcr_duplicates_downsample/{run}/{cell}.tsv"
#     threads:
#         8
#     shell:
#         """
#         sstools StatDuplicate -p {threads} -d 20 {input.bam} {output.tsv}
#         """
