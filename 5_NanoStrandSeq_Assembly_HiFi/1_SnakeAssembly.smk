#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = "results/assembly"

rule all:
    input:
        # outdir + "/flye",
        outdir + "/wtdbg2",
        # outdir + "/canu",


rule flye:
    input:
        fqs = hifi_fastqs
    output:
        out = directory(outdir + "/flye")
    log:
        outdir + "/flye.log"
    threads:
        48
    shell:
        """
        flye --pacbio-hifi  {input.fqs} --out-dir {output.out} --threads {threads} --asm-coverage 40 --genome-size 3.1g &> {log}
        """

rule wtdbg2:
    input:
        fqs = hifi_fastqs
    output:
        out = directory(outdir + "/wtdbg2")
    log:
        outdir + "/wtdbg2.log"
    threads:
        48
    shell:
        """(
        mkdir -p {output.out}
        wtdbg2 -x ccs -g 3.1g -t {threads} -o {output.out}/asm -i {input.fqs}
        wtpoa-cns -t {threads} -i {output.out}/asm.ctg.lay.gz -o {output.out}/asm.ctg.fa
        minimap2 -t {threads} -ax map-pb -r2k {output.out}/asm.ctg.fa {input.fqs} \
            | samtools sort -@ {threads} -T {output.out}/asm.ctg.bam_tmp > {output.out}/asm.ctg.bam
        samtools index -@ {threads} {output.out}/asm.ctg.bam
        samtools view -@ {threads} -F 0x900 {output.out}/asm.ctg.bam \
            | wtpoa-cns -t {threads} -d {output.out}/asm.ctg.fa -i - -o {output.out}/asm.cns.fa ) &> {log}
        """

rule canu:
    input:
        fqs = hifi_fastqs,
        out = directory(outdir + "/canu")
    log:
        outdir + "/canu.log"
    threads:
        48
    shell:
        """
        canu useGrid=false -p asm -d {output.out} genomeSize=3.1g -pacbio-hifi {input.fqs} &> {log}
        """