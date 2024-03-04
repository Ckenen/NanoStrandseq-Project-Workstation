#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = assembly_dir + "/sv"
hps = ["hp1", "hp2"]

rule all:
    input:
        outdir + "/cutesv.vcf.gz",
        outdir + "/cutesv.filtered.vcf.gz",
        outdir + "/quantify_merged.tsv",
        outdir + "/quantify_merged_lite.tsv",
        expand(outdir + "/quantify/{chrom}.{hp}.tsv.gz", chrom=chroms, hp=hps),
        expand(outdir + "/quantify_lite/{hp}.tsv", hp=hps),
        outdir + "/quantify_lite.tsv"

rule cutesv:
    input:
        fasta = lambda wildcards: GENOMES[species]["GENOME_FASTA"],
        bam = assembly_dir + "/prepare/all_cells.all_chroms.bam",
    output:
        wd = temp(directory(outdir + "/cutesv.wd")),
        vcf = temp(outdir + "/cutesv.vcf"),
        vcf2 = outdir + "/cutesv.vcf.gz"
    log:
        outdir + "/cutesv.log"
    threads:
        24 # threads
    shell:
        """(
        mkdir -p {output.wd}
        cuteSV --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads {threads} \
            --sample cuteSV \
            --retain_work_dir \
            --report_readid \
            --min_support 1 \
            --min_size 50 \
            --genotype \
            {input.bam} {input.fasta} {output.vcf} {output.wd} 
        bgzip -c {output.vcf} > {output.vcf2}
        tabix -p vcf {output.vcf2} ) &> {log}
        """

rule filter_vcf:
    input:
        vcf = rules.cutesv.output.vcf2
    output:
        vcf = outdir + "/cutesv.filtered.vcf.gz"
    shell:
        """
        zcat {input.vcf} | ../6_nss-pseudobulk-analysis/scripts/filter_cutesv_vcf.py | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
	
# pseudobulk

rule quantify_all:
    input:
        vcf = rules.filter_vcf.output.vcf,
        bam = assembly_dir + "/prepare/all_cells.all_chroms.bam"
    output:
        tsv1 = outdir + "/quantify_merged.tsv",
        tsv2 = outdir + "/quantify_merged_lite.tsv"
    threads:
        24
    shell:
        """
        sstools QuantifySV -t {threads} {input.vcf} {input.bam} {output.tsv1}
        awk -v OFS='\\t' '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}' {output.tsv1} > {output.tsv2}
        """

# haplotype-resolved

rule quantify_sv:
    input:
        vcf = rules.filter_vcf.output.vcf,
        bam = assembly_dir + "/round2/bams/{chrom}.{hp}.bam"
    output:
        vcf = temp(outdir + "/quantify/{chrom}.{hp}.vcf.gz"),
        tmp = temp(outdir + "/quantify/{chrom}.{hp}.tsv"),
        txt = outdir + "/quantify/{chrom}.{hp}.tsv.gz"
    shell:
        """
        bcftools view {input.vcf} {wildcards.chrom} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        ../6_nss-pseudobulk-analysis/scripts/quantify_sv.py -o {output.tmp} {output.vcf} {input.bam}
        pigz -p {threads} -c {output.tmp} > {output.txt}
        rm {output.vcf}.tbi
        """

rule simplify_quant_tsv:
    input:
        txt_list = [outdir + "/quantify/%s.{hp}.tsv.gz" % c for c in chroms]
    output:
        txt = outdir + "/quantify_lite/{hp}.tsv"
    run:
        import gzip
        with open(output.txt, "w+") as fw:
            for i, path in enumerate(input.txt_list):
                with gzip.open(path, "rt") as f:
                    for j, line in enumerate(f):
                        if j == 0 and i != 0:
                            continue
                        row = line.strip("\n").split("\t")[:13]
                        fw.write("\t".join(row) + "\n")

rule merge_haplotype_quant_tsv:
    input:
        txt1 = outdir + "/quantify_lite/hp1.tsv",
        txt2 = outdir + "/quantify_lite/hp2.tsv"
    output:
        txt = outdir + "/quantify_lite.tsv"
    run:
        import pandas as pd
        d1 = pd.read_csv(input.txt1, sep="\t")
        d2 = pd.read_csv(input.txt2, sep="\t")
        d1.index = d1["Name"]
        d2.index = d2["Name"]
        d = d1.merge(d2, left_index=True, right_index=True, suffixes=["_HP1", "_HP2"])
        d.to_csv(output.txt, sep="\t", index=False)
