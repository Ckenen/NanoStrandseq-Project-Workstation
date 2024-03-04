#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"

names = list(filter(lambda x: "RmDupFlag" not in x, names))

rule all:
    input:
        expand(outdir + "/bw/{name}.bw", name=names),
        expand(outdir + "/base_depth/{name}.tsv", name=names),
        # expand(outdir + "/snv_depth/{name}.json", name=names),

rule bam_to_bw:
    input:
        bam = outdir + "/bams/{name}.bam",
        txt = SIZES
    output:
        td = temp(directory(outdir + "/bw/{name}.TRACKS_TMP_TD")),
        bg = temp(outdir + "/bw/{name}.bedGraph"),
        bw = outdir + "/bw/{name}.bw"
    threads:
        12
    shell:
        """
        mkdir -p {output.td}
        bedtools genomecov -bg -ignoreD -ibam {input.bam} \
            | sort -T {output.td} --parallel {threads} -k1,1 -k2,2n > {output.bg}
        bedGraphToBigWig {output.bg} {input.txt} {output.bw}
        """

rule stat_base_depth:
    input:
        bw = outdir + "/bw/{name}.bw"
    output:
        txt = outdir + "/base_depth/{name}.tsv"
    threads:
        12
    shell:
        """
        ./scripts/stat_base_depth.py {input.bw} {threads} {output.txt}
        """

# rule stat_snv_depth:
#     input:
#         vcf = BENCHMARK_VCF,
#         bed = BENCHMARK_BED,
#         bw = outdir + "/bw/{name}.bw"
#     output:
#         txt = outdir + "/snv_depth/{name}.json"
#     threads:
#         threads
#     shell:
#         """
#         ./scripts/stat_snv_depth.py {input.vcf} {input.bed} {input.bw} {threads} {output.txt}
#         """
