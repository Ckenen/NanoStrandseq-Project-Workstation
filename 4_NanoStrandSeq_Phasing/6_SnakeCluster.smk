#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
OUTDIR = ROOT_DIR + "/clusters"

rule all:
    input:
        OUTDIR + "/anchors.bed.gz",
        expand(OUTDIR + "/assigned/{cell}.bed.gz", cell=CELLS),
        OUTDIR + "/clustered",


rule get_anchors:
    input:
        vcf = ROOT_DIR + "/snv/concat/nanocaller.vcf.gz"
    output:
        bed = OUTDIR + "/anchors.bed.gz"
    shell:
        """
        zcat {input.vcf} | grep -v '#' | grep -v '1/1' | awk 'length($4)==1 && length($5)==1' \
            | awk -v FS='\\t' -v OFS='\\t' '{{print $1,$2-1,$2,$4"/"$5}}' \
            | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """

rule assign_anchor:
    input:
        bam = ROOT_DIR + "/prepare/bams/{cell}.bam",
        bed = rules.get_anchors.output.bed,
        txt = ROOT_DIR + "/wc/{cell}/duplicate_set_names.json"
    output:
        bed = OUTDIR + "/assigned/{cell}.bed.gz"
    log:
        OUTDIR + "/assigned/{cell}.log"
    shell:
        """(
        ./scripts/strandtools/assign_anchor.py {input.bam} {input.bed} {input.txt} \
            | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
        tabix -p bed {output.bed} ) &> {log}
        """

rule clustering:
    input:
        beds = expand(OUTDIR + "/assigned/{cell}.bed.gz", cell=CELLS),
    output:
        out = directory(OUTDIR + "/clustered")
    log:
        OUTDIR + "/clustered.log"
    shell:
        """
        ./scripts/strandtools/cluster_cells.py {CONF} {output.out} &> {log}
        """