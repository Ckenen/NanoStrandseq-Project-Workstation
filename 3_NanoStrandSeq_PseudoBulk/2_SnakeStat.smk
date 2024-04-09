#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
NAMES = list(filter(lambda x: "RmDupFlag" not in x, NAMES))

rule all:
    input:
        expand(OUTDIR + "/bw/{name}.bw", name=NAMES),
        expand(OUTDIR + "/base_depth/{name}.tsv", name=NAMES),
        # expand(OUTDIR + "/snv_depth/{name}.json", name=NAMES),

rule bam_to_bw:
    input:
        bam = OUTDIR + "/bams/{name}.bam",
        txt = SIZES
    output:
        td = temp(directory(OUTDIR + "/bw/{name}.TRACKS_TMP_TD")),
        bg = temp(OUTDIR + "/bw/{name}.bedGraph"),
        bw = OUTDIR + "/bw/{name}.bw"
    threads:
        THREADS
    shell:
        """
        mkdir -p {output.td}
        bedtools genomecov -bg -ignoreD -ibam {input.bam} \
            | sort -T {output.td} --parallel {threads} -k1,1 -k2,2n > {output.bg}
        bedGraphToBigWig {output.bg} {input.txt} {output.bw}
        """

rule stat_base_depth:
    input:
        bw = rules.bam_to_bw.output.bw
    output:
        txt = OUTDIR + "/base_depth/{name}.tsv"
    threads:
        THREADS
    shell:
        """
        ./scripts/stat_base_depth.py {input.bw} {threads} {output.txt}
        """

rule stat_snv_depth:
    input:
        vcf = SNP_VCF,
        bed = SNP_BED,
        bw = rules.bam_to_bw.output.bw
    output:
        txt = OUTDIR + "/snv_depth/{name}.json"
    threads:
        THREADS
    shell:
        """
        ./scripts/stat_snv_depth.py {input.vcf} {input.bed} {input.bw} {threads} {output.txt}
        """
