#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"

# Map contigs to reference

outdir = "results/contig2genome"

rule all:
    input:
        outdir + "/minimap2/genome.map-hifi.mmi",
        outdir + "/minimap2/mapped.bam",
        outdir + "/mapped_regions.tsv",

rule build_index:
    input:
        fa = genome_fasta
    output:
        mmi = outdir + "/minimap2/genome.map-hifi.mmi"
    log:
        outdir + "/minimap2/genome.map-hifi.log"
    threads:
        24
    shell:
        """
        minimap2 -t {threads} -x map-hifi -d {output.mmi} {input.fa} &> {log}
        """

rule minimap2:
    input:
        mmi = rules.build_index.output.mmi,
        fa = contig_fasta
    output:
        bam = outdir + "/minimap2/mapped.bam"
    log:
        outdir + "/minimap2/mapped.log"
    threads:
        24
    shell:
        """(
        minimap2 -ax map-hifi --MD -t {threads} {input.mmi} {input.fa} \
            | samtools view -@ {threads} -u -F 4 - \
            | samtools sort -@ {threads} -T {output.bam}_TMP -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule get_mapped_regions:
    input:
        bam = rules.minimap2.output.bam
    output:
        txt = outdir + "/mapped_regions.tsv"
    shell:
        """
        ./scripts/get_mapped_regions.py {input.bam} {output.txt}
        """
