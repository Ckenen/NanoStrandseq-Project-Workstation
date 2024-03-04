#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = assembly_dir + "/snv"
# clair2 cost many time.
callers = ["clair2", "longshot", "nanocaller"][-1:]
# vcfeval engine need rtg, which not exists.
engines = ["xcmp", "vcfeval"][:1]

rule all:
    input:
        # expand(outdir + "/{caller}/{chrom}.vcf.gz", caller=callers, chrom=chroms),
        expand(outdir + "/concat/{caller}.vcf.gz", caller=callers),
        #expand(outdir + "/compare/{caller}.txt", caller=callers),
        #expand(outdir + "/benchmark/{engine}/{caller}", engine=engines, caller=callers),
        # outdir + "/snvs_merged.vcf.gz",
        # outdir + "/snvs_merged_common.vcf.gz",
        # outdir + "/snvs_het.bed.gz",
        
rule clair2:
    input:
        bam = assembly_dir + "/prepare/all_cells.all_chroms.bam",
        fasta = lambda wildcards: GENOMES[species]["GENOME_FASTA"]
    output:
        tmp = temp(outdir + "/clair2/{chrom}.vcf"),
        vcf = outdir + "/clair2/{chrom}.vcf.gz"
    log:
        outdir + "/clair2/{chrom}.log"
    threads:
        12
    shell:
        """(
        set +u; source activate clair
        clair.py callVarBam \
            --delay 0 \
            --chkpnt_fn /home/chenzonggui/software/princess/bin/modules/ont/model \
            --ref_fn {input.fasta} \
            --bam_fn {input.bam} \
            --ctgName {wildcards.chrom} \
            --sampleName clair2 \
            --threads {threads} \
            --call_fn {output.tmp}
        conda deactivate
        bgzip -c {output.tmp} > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule nanocaller:
    input:
        bam = assembly_dir + "/prepare/all_cells.all_chroms.bam",
        fasta = lambda wildcards: GENOMES[species]["GENOME_FASTA"]
    output:
        out = directory(outdir + "/nanocaller/{chrom}"),
        vcf = outdir + "/nanocaller/{chrom}.vcf.gz"
    log:
        outdir + "/nanocaller/{chrom}.log"
    threads:
        8
    shell:
        """(
        set +u; source activate NanoCaller
        python /home/chenzonggui/software/NanoCaller-master/scripts/NanoCaller.py \
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
        conda deactivate
        cp {output.out}/variant_calls.snps.vcf.gz {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule longshot:
    input:
        bam = assembly_dir + "/prepare/all_cells.all_chroms.bam",
        fasta = lambda wildcards: GENOMES[species]["GENOME_FASTA"]
    output:
        tmp = temp(outdir + "/longshot/{chrom}.vcf"),
        vcf = outdir + "/longshot/{chrom}.vcf.gz"
    log:
        outdir + "/longshot/{chrom}.log"
    shell:
        """(
        set +u; source activate longshot
        longshot -A \
            -r {wildcards.chrom} \
            --sample_id longshot \
            --bam {input.bam} \
            --ref {input.fasta} \
            --out {output.tmp}
        conda deactivate 
        bgzip -c {output.tmp} > {output.vcf} 
        tabix -p vcf {output.vcf} ) &> {log}
        """

def get_chrom_vcfs(wildcards):
    return [outdir + "/%s/%s.vcf.gz" % (wildcards.caller, c) for c in chroms]
        
rule concat_chrom_vcfs:
    input:
        vcfs = lambda wildcards: get_chrom_vcfs(wildcards),
    output:
        vcf = outdir + "/concat/{caller}.vcf.gz"
    log:
        outdir + "/concat/{caller}.log"
    shell:
        """(
        bcftools concat -a {input.vcfs} | bcftools sort | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule merge_caller_vcfs:
    input:
        vcf1 = outdir + "/concat/clair2.vcf.gz",
        vcf2 = outdir + "/concat/nanocaller.vcf.gz",
        vcf3 = outdir + "/concat/longshot.vcf.gz"
    output:
        vcf = outdir + "/snvs_merged.vcf.gz"
    log:
        outdir + "/snvs_merged.log"
    shell:
        """(
        bcftools merge {input.vcf1} {input.vcf2} {input.vcf3} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

# benchmark

rule compare_het_snv:
    input:
        vcf1 = lambda wildcards: BENCHMARKS[cellline]["VCF"],
        vcf2 = outdir + "/concat/{caller}.vcf.gz",
        bed = lambda wildcards: BENCHMARKS[cellline]["BED"]
    output:
        txt =  outdir + "/compare/{caller}.txt"
    shell:
        """
        ./scripts/strandtools/compare_het_snp.py {input.vcf1} {input.vcf2} {input.bed} {output.txt}
        """

rule benchmark:
    input:
        vcf1 = lambda wildcards: BENCHMARKS[cellline]["VCF"],
        vcf2 = rules.concat_chrom_vcfs.output.vcf,
        bed = lambda wildcards: BENCHMARKS[cellline]["BED"],
        fasta = lambda wildcards: GENOMES[species]["GENOME_FASTA"]
    output:
        out = directory(outdir + "/benchmark/{engine}/{caller}")
    log:
        outdir + "/benchmark/{engine}/{caller}.log"
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
        vcf = outdir + "/snvs_merged_common.vcf.gz",
        tbi = outdir + "/snvs_merged_common.vcf.gz.tbi"
    shell:
        """
        zcat {input.vcf} | grep -v '\./\.' | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule get_anchors:
    input:
        vcf = rules.merge_caller_vcfs.output.vcf
    output:
        bed = outdir + "/snvs_het.bed.gz",
        tbi = outdir + "/snvs_het.bed.gz.tbi",
    shell:
        """
        ./scripts/assembly/get_anchors.py {input.vcf} | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """


