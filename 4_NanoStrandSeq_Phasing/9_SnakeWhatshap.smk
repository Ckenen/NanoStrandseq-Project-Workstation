#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = assembly_dir

rule all:
    input:
        outdir + "/whatshap/unphase.vcf.gz",
        expand(outdir + "/whatshap/unphase/{chrom}.vcf.gz", chrom=chroms),
        #expand(outdir + "/whatshap/phased_pbccs/{chrom}.vcf.gz", chrom=chroms),
        expand(outdir + "/whatshap/phased_nss/{chrom}.vcf.gz", chrom=chroms),
        outdir + "/whatshap/phased_nss.vcf.gz",
        outdir + "/whatshap/phased_nss.benchmark_phasing.json",
        outdir + "/whatshap/phased_nss.compare.tsv"

rule unphase:
    input:
        vcf = lambda wildcards: BENCHMARKS[cellline]["VCF"]
    output:
        vcf = temp(outdir + "/whatshap/unphase.vcf"),
        vcf_gz = outdir + "/whatshap/unphase.vcf.gz"
    shell:
        """
        whatshap unphase {input.vcf} > {output.vcf}
        (
            cat {output.vcf} | grep '##'
            echo -e "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tNSS"
            cat {output.vcf} | grep -v '#'
        ) | bgzip -c > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz}
        """

rule unphase_chrom_vcf:
    input:
        vcf = outdir + "/whatshap/unphase.vcf.gz"
    output:
        vcf = outdir + "/whatshap/unphase/{chrom}.vcf.gz"
    shell:
        """
        bcftools view {input.vcf} {wildcards.chrom} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# rule phase_pbccs:
#     input:
#         fa = lambda wildcards: GENOMES[species]["GENOME_FASTA"],
#         vcf1 = outdir + "/whatshap/unphase/{chrom}.vcf.gz",
#         vcf2 = outdir + "/round2/snvs.vcf.gz",
#         bam = "/date/chenzonggui/baixiuzhen/GIAB/HG001_GRCh38_NISTv4.2.1/PacBio_SequelII_CCS_11kb/HG001_GRCh38.haplotag.RTG.trio.bam"
#     output:
#         vcf = outdir + "/whatshap/phased_pbccs/{chrom}.vcf.gz"
#     log:
#         outdir + "/whatshap/phased_pbccs/{chrom}.log"
#     shell:
#         """(
#         whatshap phase --ignore-read-groups -r {input.fa} --chromosome {wildcards.chrom} \
#             {input.vcf1} {input.vcf2} {input.bam} | bgzip -c > {output.vcf}
#         tabix -p vcf {output.vcf} ) &> {log}
#         """

rule phase_nss:
    input:
        fa = lambda wildcards: GENOMES[species]["GENOME_FASTA"],
        vcf1 = outdir + "/whatshap/unphase/{chrom}.vcf.gz",
        vcf2 = outdir + "/round2/snvs.vcf.gz",
        bam = outdir + "/prepare/all_cells.all_chroms.bam"
    output:
        vcf = outdir + "/whatshap/phased_nss/{chrom}.vcf.gz"
    log:
        outdir + "/whatshap/phased_nss/{chrom}.log"
    shell:
        """(
        whatshap phase --ignore-read-groups -r {input.fa} --chromosome {wildcards.chrom} \
            {input.vcf1} {input.vcf2} {input.bam} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule concat_chrom_vcfs:
    input:
        vcfs = expand(outdir + "/whatshap/phased_nss/{chrom}.vcf.gz", chrom=chroms)
    output:
        vcf = outdir + "/whatshap/phased_nss.vcf.gz"
    log:
        outdir + "/whatshap/phased_nss.log"
    shell:
        """(
        bcftools concat -a {input.vcfs} | bcftools sort | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule benchmark_phasing:
    input:
        vcf1 = lambda wildcards: BENCHMARKS[cellline]["VCF"],
        vcf2 = outdir + "/whatshap/phased_nss.vcf.gz",
        bed = lambda wildcards: BENCHMARKS[cellline]["BED"]
    output:
        txt = outdir + "/whatshap/phased_nss.benchmark_phasing.json"
    log:
        outdir + "/whatshap/phased_nss.benchmark_phasing.log"
    threads:
        12
    shell:
        """
        ./scripts/strandtools/benchmark_phasing.py {input} {threads} {output} &> {log}
        """

rule compare:
    input:
        vcf1 = lambda wildcards: BENCHMARKS[cellline]["VCF"],
        vcf2 = outdir + "/whatshap/phased_nss.vcf.gz"
    output:
        txt = outdir + "/whatshap/phased_nss.compare.tsv"
    log:
        outdir + "/whatshap/phased_nss.compare.log"
    shell:
        """
        whatshap compare --only-snvs --ignore-sample-name --tsv-pairwise \
            {output.txt} {input.vcf1} {input.vcf2} &> {log}
        """