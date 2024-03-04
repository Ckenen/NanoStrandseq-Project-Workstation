#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/final"
outdir = "results/stat"

rule all:
    input:
        expand(outdir + "/lengths/{run_cell}.tsv.gz", run_cell=run_cells),
        expand(outdir + "/gc_content/{run_cell}.tsv.gz", run_cell=run_cells),
        expand(outdir + "/gc_density/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/background/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/spikiness/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/depth/{run_cell}.tsv", run_cell=run_cells),
        expand(outdir + "/coverage/{run_cell}.tsv", run_cell=run_cells),

rule stat_length:
    input:
        bam = indir + "/{run}/{cell}.bam"
    output:
        tmp = temp(outdir + "/lengths/{run}/{cell}.tsv"),
        tsv = outdir + "/lengths/{run}/{cell}.tsv.gz",
        tsv2 = outdir + "/lengths/{run}/{cell}_summary.tsv"
    log:
        outdir + "/lengths/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatLength -t {threads} -l SE -s {output.tsv2} {input.bam} {output.tmp} &> {log}
        gzip -c {output.tmp} > {output.tsv}
        """

rule stat_gc_content:
    input:
        bam = indir + "/{run}/{cell}.bam",
        fa = config["fasta"]
    output:
        tmp = temp(outdir + "/gc_content/{run}/{cell}.tsv"),
        tsv = outdir + "/gc_content/{run}/{cell}.tsv.gz",
        tsv2 = outdir + "/gc_content/{run}/{cell}_summary.tsv"
    log:
        outdir + "/gc_content/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatGC -t {threads} -s {output.tsv2} {input.bam} {input.fa} {output.tmp} &> {log}
        gzip -c {output.tmp} > {output.tsv}
        """

rule prepare_genome_gc:
    input:
        fa = config["fasta"]
    output:
        pkl = outdir + "/gc_density_pkl/%s.pkl" % config["build"]
    log:
        outdir + "/gc_density_pkl/%s.log" % config["build"]
    threads:
        12
    shell:
        """
        sstools StatGCDensity --prepare -t {threads} -k {output.pkl} {input.fa} &> {log}
        """

rule stat_gc_density:
    input:
        bam = indir + "/{run}/{cell}.bam",
        fa = config["fasta"],
        pkl = outdir + "/gc_density_pkl/%s.pkl" % config["build"]
    output:
        tsv = outdir + "/gc_density/{run}/{cell}.tsv"
    log:
        outdir + "/gc_density/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatGCDensity -t {threads} -k {input.pkl} {input.bam} {input.fa} {output.tsv} &> {log}
        """

rule stat_background:
    input:
        bam = indir + "/{run}/{cell}.bam"
    output:
        tsv = outdir + "/background/{run}/{cell}.tsv",
        tsv2 = outdir + "/background/{run}/{cell}_summary.tsv"
    log:
        outdir + "/background/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatBackground -t {threads} -s {output.tsv2} {input.bam} {output.tsv} &> {log}
        """

rule stat_spikiness:
    input:
        bam = indir + "/{run}/{cell}.bam"
    output:
        tsv = outdir + "/spikiness/{run}/{cell}.tsv"
    log:
        outdir + "/spikiness/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatSpikiness -t {threads} {input.bam} {output.tsv} &> {log}
        """

rule stat_depth:
    input:
        bam = indir + "/{run}/{cell}.bam"
    output:
        tsv1 = outdir + "/depth/{run}/{cell}.tsv",
        tsv2 = outdir + "/depth/{run}/{cell}_rmdup.tsv"
    log:
        outdir + "/depth/{run}/{cell}.log"
    threads:
        4
    shell:
        """(
        sstools StatDepth -t {threads} {input.bam} {output.tsv1}
        sstools StatDepth -t {threads} -r {input.bam} {output.tsv2} ) &> {log}
        """

rule stat_coverage:
    input:
        bam = indir + "/{run}/{cell}.bam"
    output:
        tsv = outdir + "/coverage/{run}/{cell}.tsv"
    log:
        outdir + "/coverage/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools StatCoverage -t {threads} {input.bam} {output.tsv} &> {log}
        """

# rule stat_genomic_coverage_downsample:
#     input:
#         bam = indir + "/{run}/{cell}.bam"
#     output:
#         tsv = outdir + "/genomic_coverage_downsample/{run}/{cell}.tsv"
#     threads:
#         8
#     shell:
#         """
#         sstools StatGenomicCov -p {threads} {input.bam} {output.tsv}
#         """


# rule stat_pcr_duplicates_downsample:
#     input:
#         bam = indir + "/{run}/{cell}.bam",
#     output:
#         tsv = outdir + "/pcr_duplicates_downsample/{run}/{cell}.tsv"
#     threads:
#         8
#     shell:
#         """
#         sstools StatDuplicate -p {threads} -d 20 {input.bam} {output.tsv}
#         """
