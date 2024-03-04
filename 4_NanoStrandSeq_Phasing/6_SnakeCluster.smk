#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = assembly_dir + "/clusters"

rule all:
    input:
        outdir + "/anchors.bed.gz",
        expand(outdir + "/assigned/{cell}.bed.gz", cell=cells),
        outdir + "/clustered",


rule get_anchors:
    input:
        vcf = assembly_dir + "/snv/concat/nanocaller.vcf.gz"
    output:
        bed = outdir + "/anchors.bed.gz"
    shell:
        """
        zcat {input.vcf} | grep -v '#' | grep -v '1/1' | awk 'length($4)==1 && length($5)==1' \
            | awk -v FS='\\t' -v OFS='\\t' '{{print $1,$2-1,$2,$4"/"$5}}' \
            | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """

rule assign_anchor:
    input:
        bam = assembly_dir + "/prepare/bams/{cell}.bam",
        bed = rules.get_anchors.output.bed,
        txt = assembly_dir + "/wc/{cell}/duplicate_set_names.json"
    output:
        bed = outdir + "/assigned/{cell}.bed.gz"
    log:
        outdir + "/assigned/{cell}.log"
    shell:
        """(
        ./scripts/strandtools/assign_anchor.py {input.bam} {input.bed} {input.txt} \
            | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
        tabix -p bed {output.bed} ) &> {log}
        """

rule clustering:
    input:
        beds = expand(outdir + "/assigned/{cell}.bed.gz", cell=cells),
    output:
        out = directory(outdir + "/clustered")
    log:
        outdir + "/clustered.log"
    shell:
        """
        ./scripts/strandtools/cluster_cells.py {assembly_conf} {output.out} &> {log}
        """