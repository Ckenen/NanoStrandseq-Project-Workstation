#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
OUTDIR = ROOT_DIR + "/snv"
# clair2 cost many time.
callers = ["clair2", "longshot", "nanocaller"][-1:]
# vcfeval engine need rtg, which not exists.
engines = ["xcmp", "vcfeval"][:1]

rule all:
    input:
        # expand(OUTDIR + "/{caller}/{chrom}.vcf.gz", caller=callers, chrom=CHROMS),
        expand(OUTDIR + "/concat/{caller}.vcf.gz", caller=callers),
        #expand(OUTDIR + "/compare/{caller}.txt", caller=callers),
        #expand(OUTDIR + "/benchmark/{engine}/{caller}", engine=engines, caller=callers),
        # OUTDIR + "/snvs_merged.vcf.gz",
        # OUTDIR + "/snvs_merged_common.vcf.gz",
        # OUTDIR + "/snvs_het.bed.gz",
        
rule clair2:
    input:
        bam = ROOT_DIR + "/prepare/all_cells.all_chroms.bam",
        fasta = GENOME_FASTA
    output:
        tmp = temp(OUTDIR + "/clair2/{chrom}.vcf"),
        vcf = OUTDIR + "/clair2/{chrom}.vcf.gz"
    log:
        OUTDIR + "/clair2/{chrom}.log"
    conda:
        "clair"
    threads:
        THREADS
    shell:
        """(
        clair.py callVarBam \
            --delay 0 \
            --chkpnt_fn ${{HOME}}/software/princess/bin/modules/ont/model \
            --ref_fn {input.fasta} \
            --bam_fn {input.bam} \
            --ctgName {wildcards.chrom} \
            --sampleName clair2 \
            --threads {threads} \
            --call_fn {output.tmp}
        bgzip -c {output.tmp} > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule nanocaller:
    input:
        bam = ROOT_DIR + "/prepare/all_cells.all_chroms.bam",
        fasta = GENOME_FASTA
    output:
        out = directory(OUTDIR + "/nanocaller/{chrom}"),
        vcf = OUTDIR + "/nanocaller/{chrom}.vcf.gz"
    log:
        OUTDIR + "/nanocaller/{chrom}.log"
    conda:
        "NanoCaller"
    threads:
        THREADS
    shell:
        """(
        python ${{HOME}}/software/NanoCaller-master/scripts/NanoCaller.py \
            --mode snps_unphased \
            -bam {input.bam} \
            -p ont \
            -o {output.out} \
            -chrom {wildcards.chrom} \
            -ref {input.fasta} \
            -cpu {threads} \
            --snp_model ONT-HG001 \
            --indel_model ONT-HG001 \
            -sample nanocaller
        cp {output.out}/variant_calls.snps.vcf.gz {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule longshot:
    input:
        bam = ROOT_DIR + "/prepare/all_cells.all_chroms.bam",
        fasta = lambda wildcards: GENOME_FASTA
    output:
        tmp = temp(OUTDIR + "/longshot/{chrom}.vcf"),
        vcf = OUTDIR + "/longshot/{chrom}.vcf.gz"
    log:
        OUTDIR + "/longshot/{chrom}.log"
    conda:
        "longshot"
    shell:
        """(
        longshot -A \
            -r {wildcards.chrom} \
            --sample_id longshot \
            --bam {input.bam} \
            --ref {input.fasta} \
            --out {output.tmp}
        bgzip -c {output.tmp} > {output.vcf} 
        tabix -p vcf {output.vcf} ) &> {log}
        """

def get_chrom_vcfs(wildcards):
    return [OUTDIR + "/%s/%s.vcf.gz" % (wildcards.caller, c) for c in CHROMS]
        
rule concat_chrom_vcfs:
    input:
        vcfs = lambda wildcards: get_chrom_vcfs(wildcards),
    output:
        vcf = OUTDIR + "/concat/{caller}.vcf.gz"
    log:
        OUTDIR + "/concat/{caller}.log"
    shell:
        """(
        bcftools concat -a {input.vcfs} | bcftools sort | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule merge_caller_vcfs:
    input:
        vcf1 = OUTDIR + "/concat/clair2.vcf.gz",
        vcf2 = OUTDIR + "/concat/nanocaller.vcf.gz",
        vcf3 = OUTDIR + "/concat/longshot.vcf.gz"
    output:
        vcf = OUTDIR + "/snvs_merged.vcf.gz"
    log:
        OUTDIR + "/snvs_merged.log"
    shell:
        """(
        bcftools merge {input.vcf1} {input.vcf2} {input.vcf3} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

# benchmark

rule compare_het_snv:
    input:
        vcf1 = SNP_VCF,
        vcf2 = OUTDIR + "/concat/{caller}.vcf.gz",
        bed = SNP_BED
    output:
        txt =  OUTDIR + "/compare/{caller}.txt"
    shell:
        """
        ./scripts/strandtools/compare_het_snp.py {input.vcf1} {input.vcf2} {input.bed} {output.txt}
        """

rule benchmark:
    input:
        vcf1 = SNP_VCF,
        vcf2 = rules.concat_chrom_vcfs.output.vcf,
        bed = SNP_BED,
        fasta = GENOME_FASTA
    output:
        out = directory(OUTDIR + "/benchmark/{engine}/{caller}")
    log:
        OUTDIR + "/benchmark/{engine}/{caller}.log"
    threads:
        8
    shell:
        """(
        mkdir -p {output.out}
        set +u; source activate happy
        hap.py -r {input.fasta} \
            -o {output.out}/benchmark \
            --threads {threads} \
            --engine {wildcards.engine} \
            -T {input.bed} \
            {input.vcf1} {input.vcf2} ) &> {log}
        """

rule get_common_vcfs:
    input:
        vcf = rules.merge_caller_vcfs.output.vcf
    output:
        vcf = OUTDIR + "/snvs_merged_common.vcf.gz",
        tbi = OUTDIR + "/snvs_merged_common.vcf.gz.tbi"
    shell:
        """
        zcat {input.vcf} | grep -v '\./\.' | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule get_anchors:
    input:
        vcf = rules.merge_caller_vcfs.output.vcf
    output:
        bed = OUTDIR + "/snvs_het.bed.gz",
        tbi = OUTDIR + "/snvs_het.bed.gz.tbi",
    shell:
        """
        ./scripts/assembly/get_anchors.py {input.vcf} | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """


