#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"

# For NanoStrand-seq
import pandas as pd
select_cells = pd.read_csv("data/select_nanostrandseq_cells.20230925.tsv", sep="\t")

rule all:
    input:
        expand(outdir + "/bams/{name}.bam", name=names),
        expand(outdir + "/bams/{name}.flagstat", name=names),
        expand(outdir + "/genome_depth/{name}.json", name=names),

# Full PB-CCS, ONT-UL, and NSS

rule make_pbccs_full_bam:
    input:
        bam = "data/HG001_GRCh38.haplotag.RTG.trio.bam"
    output:
        bam = outdir + "/bams/PacBio.full.bam"
    threads:
        4
    shell:
        """
        ./scripts/filter_giab_bam.py {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule make_ontul_full_bam:
    input:
        bam = "data/NA12878-minion-ul_GRCh38.bam"
    output:
        bam = outdir + "/bams/Ultralong.full.bam"
    threads:
        4
    shell:
        """
        ./scripts/filter_giab_bam.py {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule make_nss_full_bam:
   output:
       bam = outdir + "/bams/NSS.full.bam"
   threads:
       12
   shell:
       """
       samtools merge -@ {threads} -o {output.bam} {NSS_BAM_DIR}/*.bam
       samtools index -@ {threads} {output.bam}
       """

# Downsample for NanoStrand-seq

def get_nss_bams(wildcards):
    name = "cov%s-r%s" % (wildcards.cov, wildcards.rep)
    cells = select_cells[select_cells["Name"] == name]["Cells"].values[0]
    cells = cells.split(",")
    paths = ["data/nss_bams/%s.bam" % c for c in cells]
    return paths

rule make_nss_downsample_bam:
    input:
        bams = lambda wildcards: get_nss_bams(wildcards)
    output:
        bam = outdir + "/bams/NSS.cov{cov}-r{rep}.bam"
    threads:
        12
    shell:
        """
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam}
        """

rule make_nss_downsample_rmdupflag_bam:
    input:
        bam = outdir + "/bams/NSS.cov{cov}-r{rep}.bam"
    output:
        bam = outdir + "/bams/NSS_RmDupFlag.cov{cov}-r{rep}.bam"
    threads:
        4
    shell:
        """
        ./scripts/rm_dup_flag.py {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam}
        """

# Downsample for PacBio, Ultralong

rule make_bulk_downsample_bam:
    input:
        bam = outdir + "/bams/{source}.full.bam",
        txt = outdir + "/genome_depth/{source}.full.json"
    output:
        bam = outdir + "/bams/{source,PacBio|Ultralong}.cov{cov}-r{rep}.bam"
    threads:
        4
    shell:
        """
        p=`grep Depth {input.txt} | head -n 1 | awk '{{print {wildcards.cov}/$2}}'`
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
        bam = outdir + "/bams/{name}.bam"
    output:
        txt = outdir + "/genome_depth/{name}.json"
    threads:
        8
    shell:
        """
        sstools CalGenomeDepth -t {threads} -o {output.txt} {input.bam}
        """