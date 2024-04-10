#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
OUTDIR = ROOT_DIR

rule all:
    input:
        OUTDIR + "/whatshap/unphase.vcf.gz",
        expand(OUTDIR + "/whatshap/unphase/{chrom}.vcf.gz", chrom=CHROMS),
        # expand(OUTDIR + "/whatshap/phased_pbccs/{chrom}.vcf.gz", chrom=CHROMS),
        expand(OUTDIR + "/whatshap/phased_nss/{chrom}.vcf.gz", chrom=CHROMS),
        OUTDIR + "/whatshap/phased_nss.vcf.gz",
        OUTDIR + "/whatshap/phased_nss.benchmark_phasing.json",
        OUTDIR + "/whatshap/phased_nss.compare.tsv"

rule unphase:
    input:
        vcf = SNP_VCF
    output:
        vcf = temp(OUTDIR + "/whatshap/unphase.vcf"),
        vcf_gz = OUTDIR + "/whatshap/unphase.vcf.gz"
    conda:
        "whatshap"
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
        vcf = OUTDIR + "/whatshap/unphase.vcf.gz"
    output:
        vcf = OUTDIR + "/whatshap/unphase/{chrom}.vcf.gz"
    shell:
        """
        bcftools view {input.vcf} {wildcards.chrom} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# rule phase_pbccs:
#     input:
#         fa = GENOME_FASTA,
#         vcf1 = OUTDIR + "/whatshap/unphase/{chrom}.vcf.gz",
#         vcf2 = OUTDIR + "/round2/snvs.vcf.gz",
#         bam = "/date/chenzonggui/baixiuzhen/GIAB/HG001_GRCh38_NISTv4.2.1/PacBio_SequelII_CCS_11kb/HG001_GRCh38.haplotag.RTG.trio.bam"
#     output:
#         vcf = OUTDIR + "/whatshap/phased_pbccs/{chrom}.vcf.gz"
#     log:
#         OUTDIR + "/whatshap/phased_pbccs/{chrom}.log"
#     shell:
#         """(
#         whatshap phase --ignore-read-groups -r {input.fa} --chromosome {wildcards.chrom} \
#             {input.vcf1} {input.vcf2} {input.bam} | bgzip -c > {output.vcf}
#         tabix -p vcf {output.vcf} ) &> {log}
#         """

rule phase_nss:
    input:
        fa = GENOME_FASTA,
        vcf1 = OUTDIR + "/whatshap/unphase/{chrom}.vcf.gz",
        vcf2 = OUTDIR + "/round2/snvs.vcf.gz",
        bam = OUTDIR + "/prepare/all_cells.all_chroms.bam"
    output:
        vcf = OUTDIR + "/whatshap/phased_nss/{chrom}.vcf.gz"
    log:
        OUTDIR + "/whatshap/phased_nss/{chrom}.log"
    conda:
        "whatshap"
    shell:
        """(
        whatshap phase --ignore-read-groups -r {input.fa} --chromosome {wildcards.chrom} \
            {input.vcf1} {input.vcf2} {input.bam} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule concat_chrom_vcfs:
    input:
        vcfs = expand(OUTDIR + "/whatshap/phased_nss/{chrom}.vcf.gz", chrom=CHROMS)
    output:
        vcf = OUTDIR + "/whatshap/phased_nss.vcf.gz"
    log:
        OUTDIR + "/whatshap/phased_nss.log"
    shell:
        """(
        bcftools concat -a {input.vcfs} | bcftools sort | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule benchmark_phasing:
    input:
        vcf1 = SNP_VCF,
        vcf2 = OUTDIR + "/whatshap/phased_nss.vcf.gz",
        bed = SNP_BED
    output:
        txt = OUTDIR + "/whatshap/phased_nss.benchmark_phasing.json"
    log:
        OUTDIR + "/whatshap/phased_nss.benchmark_phasing.log"
    threads:
        12
    shell:
        """
        ./scripts/strandtools/benchmark_phasing.py {input} {threads} {output} &> {log}
        """

rule compare:
    input:
        vcf1 = SNP_VCF,
        vcf2 = OUTDIR + "/whatshap/phased_nss.vcf.gz"
    output:
        txt = OUTDIR + "/whatshap/phased_nss.compare.tsv"
    log:
        OUTDIR + "/whatshap/phased_nss.compare.log"
    conda:
        "whatshap"
    shell:
        """
        whatshap compare --only-snvs --ignore-sample-name --tsv-pairwise \
            {output.txt} {input.vcf1} {input.vcf2} &> {log}
        """