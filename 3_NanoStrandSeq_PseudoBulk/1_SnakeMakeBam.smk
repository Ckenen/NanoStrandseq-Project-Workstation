#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
import glob

# For NanoStrand-seq
import pandas as pd
DAT = pd.read_csv("data/select_nanostrandseq_cells.tsv", sep="\t")

rule all:
    input:
        expand(OUTDIR + "/bams/{name}.bam", name=NAMES),
        expand(OUTDIR + "/bams/{name}.flagstat", name=NAMES),
        expand(OUTDIR + "/genome_depth/{name}.tsv", name=NAMES),

# Full PB-CCS, ONT-UL, and NSS

rule make_pbccs_full_bam:
    input:
        bam = PBCCS_BAM
    output:
        bam = OUTDIR + "/bams/PacBio.full.bam"
    threads:
        1
    shell:
        """
        ./scripts/filter_giab_bam.py {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule make_ontul_full_bam:
    input:
        bam = ONTUL_BAM
    output:
        bam = OUTDIR + "/bams/Ultralong.full.bam"
    threads:
        1
    shell:
        """
        ./scripts/filter_giab_bam.py {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule make_nss_full_bam:
    input:
        bams = list(sorted(glob.glob(NSS_BAM_DIR + "/*.bam")))
    output:
        bam = OUTDIR + "/bams/NSS.full.bam"
    threads:
        THREADS
    shell:
        """
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam}
        """

# Downsample for NanoStrand-seq

def get_nss_bams(wildcards):
    name = "cov%s-r%s" % (wildcards.cov, wildcards.rep)
    cells = DAT[DAT["Name"] == name]["Cells"].values[0]
    cells = cells.split(",")
    paths = ["data/nss_bams/%s.bam" % c for c in cells]
    return paths

rule make_nss_downsample_bam:
    input:
        bams = lambda wildcards: get_nss_bams(wildcards)
    output:
        bam = OUTDIR + "/bams/NSS.cov{cov}-r{rep}.bam"
    threads:
        THREADS
    shell:
        """
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam}
        """

rule make_nss_downsample_rmdupflag_bam:
    input:
        bam = OUTDIR + "/bams/NSS.cov{cov}-r{rep}.bam"
    output:
        bam = OUTDIR + "/bams/NSS_RmDupFlag.cov{cov}-r{rep}.bam"
    threads:
        1
    shell:
        """
        ./scripts/rm_dup_flag.py {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam}
        """

# Downsample for PacBio, Ultralong

rule make_bulk_downsample_bam:
    input:
        bam = OUTDIR + "/bams/{source}.full.bam",
        txt = OUTDIR + "/genome_depth/{source}.full.tsv"
    output:
        bam = OUTDIR + "/bams/{source,PacBio|Ultralong}.cov{cov}-r{rep}.bam"
    threads:
        1
    shell:
        """
        p=`head -n 1 {input.txt} | awk '{{print {wildcards.cov}/$4}}'`
        sstools DownsampleBam -s {wildcards.rep} -p $p {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam}
        """

# Common rules

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """

rule cal_genome_depth:
    input:
        bam = OUTDIR + "/bams/{name}.bam"
    output:
        tsv = OUTDIR + "/genome_depth/{name}.tsv"
    threads:
        THREADS
    shell:
        """
        sstools StatDepth -t {threads} {input.bam} {output.tsv}
        """