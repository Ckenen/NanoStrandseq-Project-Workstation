#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = assembly_dir
hps = ["hp1", "hp2"]
#chroms = ["chr1"]
# chroms.remove("chrX")
rs = ["round1", "round2"]
# chroms = ["chrX"]

rule all:
    input:
        # expand(outdir + "/{r}/bams/{chrom}.{hp}.bam", r=rs, chrom=chroms, hp=hps),
        # expand(outdir + "/{r}/bams/{chrom}.{hp}.flagstat", r=rs, chrom=chroms, hp=hps),
        expand(outdir + "/{r}/merged.bam", r=["round2"]),
        # expand(outdir + "/{r}/merged.{hp}.bam", r=rs, hp=hps),
        # expand(outdir + "/{r}/quantify_psudobulk_sv.{hp}_lite.tsv", r=rs, hp=hps),
        expand(outdir + "/{r}/matrix/{chrom}.{hp}", r=rs, chrom=chroms, hp=hps),
        #expand(outdir + "/{r}/benchmark/{chrom}.{hp}", r=rs, chrom=chroms, hp=hps),
        expand(outdir + "/{r}/matrix2/{chrom}", r=rs, chrom=chroms),
        expand(outdir + "/{r}/matrix2.stat/{chrom}.tsv", r=rs, chrom=chroms),
        expand(outdir + "/{r}/matrix2.stat.plot/{chrom}", r=rs, chrom=chroms),
        expand(outdir + "/{r}/matrix2.filtered/{chrom}.matrix.gz", r=rs, chrom=chroms),
        # expand(outdir + "/{r}/snvs/{chrom}.vcf.gz", r=rs, chrom=chroms),
        # expand(outdir + "/{r}/snvs/{chrom}.benchmark", r=rs, chrom=chroms),
        expand(outdir + "/{r}/snvs.vcf.gz", r=rs),
        expand(outdir + "/{r}/snvs_benchmark.json", r=rs),
        #expand(outdir + "/{r}/snvs.benchmark", r=rs),
        # expand(outdir + "/{r}/hets/{chrom}.bed", r=rs, chrom=chroms),
        expand(outdir + "/{r}/hets.all.bed.gz", r=rs),
        #expand(outdir + "/{r}/hets.benchmark", r=rs),
        #expand(outdir + "/{r}/cuteSV/{chrom}.{hp}.vcf.gz", r=rs, chrom=chroms, hp=hps),
        #expand(outdir + "/{r}/cuteSV.{hp}.vcf.gz", r=rs, hp=hps),
        #expand(outdir + "/{r}/cutesv/quantify/{chrom}.{hp}.on_{hpA}.txt.gz", r=rs, chrom=chroms, hp=hps, hpA=hps),
        #expand(outdir + "/{r}/cutesv/quantify_lite/{chrom}.{hp}.on_{hpA}.tsv", r=rs, chrom=chroms, hp=hps, hpA=hps),
        
## Round 1

rule make_round1_bam:
    input:
        txt = assembly_dir + "/clusters/clustered/cells.{chrom}.json"
    output:
        tmp1 = temp(outdir + "/round1/bams/{chrom}.hp1.unsorted.bam"),
        tmp2 = temp(outdir + "/round1/bams/{chrom}.hp2.unsorted.bam"),
        bam1 = outdir + "/round1/bams/{chrom}.hp1.bam",
        bam2 = outdir + "/round1/bams/{chrom}.hp2.bam"
    log:
        outdir + "/round1/bams/{chrom}.log"
    threads:
        4
    shell:
        """(
        ./scripts/strandtools/make_haplotype_bam.py \
            {input.txt} {assembly_dir} {wildcards.chrom} \
            {output.tmp1} {output.tmp2}
        samtools sort -@ {threads} -o {output.bam1} {output.tmp1}
        samtools sort -@ {threads} -o {output.bam2} {output.tmp2}
        samtools index -@ {threads} {output.bam1}
        samtools index -@ {threads} {output.bam2} ) &> {log}
        """

## Round 2

rule split_haplotype:
    input:
        bam = assembly_dir + "/prepare/bams/{cell}.bam",
        bed = outdir + "/round1/hets.all.bed.gz"
    output:
        out = directory(outdir + "/round2/splited_haplotype/{cell}")
    log:
        outdir + "/round2/splited_haplotype/{cell}.log"
    threads:
        4
    shell:
        """
        ./scripts/assembly/split_haplotype.v1.py {input.bam} {input.bed} {output.out} &> {log}
        """

rule make_round2_bam:
    input:
        expand(rules.split_haplotype.output.out, cell=cells),
    output:
        tmp1 = temp(outdir + "/round2/bams/{chrom}.hp1.tmp.bam"),
        tmp2 = temp(outdir + "/round2/bams/{chrom}.hp2.tmp.bam"),
        bam1 = outdir + "/round2/bams/{chrom}.hp1.bam",
        bam2 = outdir + "/round2/bams/{chrom}.hp2.bam"
    log:
        outdir + "/round2/bams/{chrom}.log"
    threads:
        4
    shell:
        """(
        ./scripts/strandtools/make_haplotype_bam_r2.py \
            {assembly_conf} {wildcards.chrom} {output.tmp1} {output.tmp2}
        samtools sort -@ {threads} -o {output.bam1} {output.tmp1}
        samtools sort -@ {threads} -o {output.bam2} {output.tmp2}
        samtools index -@ {threads} {output.bam1}
        samtools index -@ {threads} {output.bam2} ) &> {log}
        """

## Common steps

rule merge_bams:
    input:
        bams1 = [outdir + "/{round}/bams/%s.hp1.bam" % c for c in chroms],
        bams2 = [outdir + "/{round}/bams/%s.hp2.bam" % c for c in chroms]
    output:
        bam = outdir + "/{round}/merged.bam"
    threads:
        12
    shell:
        """
        (   
            samtools view -H --no-PG {input.bams1[0]}
            for f in {input.bams1}; do
                samtools view -@ {threads} $f | awk '{{print $0"\\tHP:Z:HP1"}}'
            done
            for f in {input.bams2}; do
                samtools view -@ {threads} $f | awk '{{print $0"\\tHP:Z:HP2"}}'
            done 
        ) | samtools view -@ {threads} -u - \
            | samtools sort -@ {threads} -T {output.bam}_TEMP - > {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule merge_haplotype_bam:
    input:
        bams = [outdir + "/{round}/bams/%s.{hp}.bam" % c for c in chroms]
    output:
        bam  = outdir + "/{round}/merged.{hp}.bam"
    threads:
        8
    shell:
        """
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam}
        """

rule generate_base_matrix:
    input:
        bam = outdir + "/{r}/bams/{chrom}.{hp}.bam",
        fa = lambda wildcards: GENOMES[species]["GENOME_FASTA"]
    output:
        mtx = directory(outdir + "/{r}/matrix/{chrom}.{hp}")
    log:
        outdir + "/{r}/matrix/{chrom}.{hp}.log"
    threads:
        24
    shell:
        """
        ./scripts/assembly/generate_base_matrix.py \
            {input.bam} {input.fa} {wildcards.chrom} {threads} {output.mtx} &> {log}
        """

rule merge_base_matrix:
    input:
        fa = lambda wildcards: GENOMES[species]["GENOME_FASTA"],
        bed = lambda wildcards: BENCHMARKS[cellline]["BED"],
        vcf = lambda wildcards: BENCHMARKS[cellline]["VCF"],
        mtxdir1 = outdir + "/{r}/matrix/{chrom}.hp1",
        mtxdir2 = outdir + "/{r}/matrix/{chrom}.hp2"
    output:
        directory(outdir + "/{r}/matrix2/{chrom}")
    log:
        log = outdir + "/{r}/matrix2/{chrom}.log"
    threads:
        12
    shell:
        """
        ./scripts/assembly/merge_base_matrix.py {input.fa} {input.bed} {input.vcf} \
            {input.mtxdir1} {input.mtxdir2} {wildcards.chrom} {threads} {output} &> {log}
        """

rule stat_matrix2:
    input:
        outdir + "/{r}/matrix2/{chrom}"
    output:
        tsv = outdir + "/{r}/matrix2.stat/{chrom}.tsv"
    log:
        outdir + "/{r}/matrix2.stat/{chrom}.log"
    threads:
        12
    shell:
        """
        ./scripts/assembly/stat_matrix2.py {input} {threads} {output.tsv} &> {log}
        """

rule plot_stat_matrix2:
    input:
        tsv = rules.stat_matrix2.output.tsv
    output:
        directory(outdir + "/{r}/matrix2.stat.plot/{chrom}")
    shell:
        """
        ./scripts/assembly/plot_stat_matrix2.py {input.tsv} {output}
        """

rule filter_matrix2:
    input:
        outdir + "/{r}/matrix2/{chrom}"
    output:
        txt = outdir + "/{r}/matrix2.filtered/{chrom}.matrix.gz"
    shell:
        """
        zcat {input}/*.matrix.gz \
            | awk '$6!="."&&$7!="."&&$6!="-"&&$7!="-"&&($6!=$2||$7!=$2)' \
            | sort -k1,1n | gzip -c > {output.txt}
        """

# VCF

rule make_chrom_snvs:
    input:
        mtx = outdir + "/{r}/matrix2.filtered/{chrom}.matrix.gz",
        vcf = assembly_dir + "/snv/concat/nanocaller.vcf.gz",
        bed = assembly_dir + "/inversions/inversions.bed.gz",
        txt = lambda wildcards: GENOMES[species]["GENOME_SIZE"]
    output:
        tmp = temp(outdir + "/{r}/snvs/{chrom}.vcf"),
        vcf = outdir + "/{r}/snvs/{chrom}.vcf.gz"
    log:
        outdir + "/{r}/snvs/{chrom}.log"
    shell:
        """(
        ./scripts/assembly/call_snps.py {input.mtx} {input.vcf} {input.bed} {input.txt} {wildcards.chrom} {output.tmp}
        bgzip -c {output.tmp} > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule concat_chrom_vcfs:
    input:
        vcfs = [outdir + "/{r}/snvs/%s.vcf.gz" % c for c in chroms]
    output:
        vcf = outdir + "/{r}/snvs.vcf.gz"
    log:
        outdir + "/{r}/snvs.log"
    shell:
        """(
        bcftools concat -a {input.vcfs} | bcftools sort | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule benchmark_snps:
    input:
        vcf1 = lambda wildcards: BENCHMARKS[cellline]["VCF"],
        vcf2 = outdir + "/{r}/snvs.vcf.gz",
        bed = lambda wildcards: BENCHMARKS[cellline]["BED"]
    output:
        txt = outdir + "/{r}/snvs_benchmark.json"
    log:
         outdir + "/{r}/snvs_benchmark.log"
    threads:
        12
    shell:
        """
        sstools BenchmarkSNP -b {input.bed} -t {threads} \
            -o {output.txt} {input.vcf1} {input.vcf2} &> {log}
        """

# hetSNPs

rule get_hets:
    input:
        bed = assembly_dir + "/snv/concat/nanocaller.vcf.gz",
        mtxdir = outdir + "/{r}/matrix2/{chrom}"
    output:
        bed = outdir + "/{r}/hets/{chrom}.bed",
        bed_gz = outdir + "/{r}/hets/{chrom}.bed.gz"
    log:
        outdir + "/{r}/hets/{chrom}.log"
    threads:
        12
    shell:
        """(
        ./scripts/assembly/get_hets.p.py {input.bed} {input.mtxdir} {wildcards.chrom} {threads} {output.bed} 
        bgzip -c {output.bed} > {output.bed_gz}
        tabix -p bed {output.bed_gz} ) &> {log}
        """

rule merge_hets:
    input:
        beds = [outdir + "/{r}/hets/%s.bed" % c for c in chroms]
    output:
        bed = outdir + "/{r}/hets.all.bed.gz"
    shell:
        """
        cat {input.beds} | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """

# SV

rule cutesv:
    input:
        bam = outdir + "/{r}/bams/{chrom}.{hp}.bam",
        fa = lambda wildcards: GENOMES[species]["GENOME_FASTA"]
    output:
        wd = temp(directory(outdir + "/{r}/cuteSV/{chrom}.{hp}.wd")),
        vcf = temp(outdir + "/{r}/cuteSV/{chrom}.{hp}.vcf"),
        vcf2 = outdir + "/{r}/cuteSV/{chrom}.{hp}.vcf.gz"
    log:
        outdir + "/{r}/cuteSV/{chrom}.{hp}.log"
    threads:
        8
    shell:
        """(
        mkdir -p {output.wd}
        cuteSV --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads {threads} \
            --sample HG001 \
            --retain_work_dir \
            --report_readid \
            --min_support 1 \
            --min_size 50 \
            {input.bam} {input.fa} \
            {output.vcf} {output.wd} 
        bgzip -c {output.vcf} > {output.vcf2}
        tabix -p vcf {output.vcf2} ) &> {log}
        """

rule concat_cutesv_vcfs:
    input:
        vcfs = [outdir + "/{r}/cuteSV/%s.{hp}.vcf.gz" % c for c in chroms]
    output:
        vcf = outdir + "/{r}/cuteSV.{hp}.vcf.gz"
    shell:
        """
        bcftools concat -a {input.vcfs} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule quantify_sv:
    input:
        vcf = outdir + "/{r}/cuteSV/{chrom}.{hp}.vcf.gz",
        bam = outdir + "/{r}/bams/{chrom}.{hpA}.bam"
    output:
        tmp = temp(outdir + "/{r}/cutesv/quantify/{chrom}.{hp}.on_{hpA}.txt"),
        txt = outdir + "/{r}/cutesv/quantify/{chrom}.{hp}.on_{hpA}.txt.gz"
    threads:
        threads
    shell:
        """
        ../6_SNV_SV_Comparison/scripts/quantify_sv.p.py {input.vcf} {threads} {input.bam} {output.tmp}
        awk 'NR==1||$1=="{wildcards.chrom}"' {output.tmp} | pigz -p {threads} -c > {output.txt}
        """

rule lite_quantify:
    input:
        txt = outdir + "/{r}/cutesv/quantify/{source}.txt.gz"
    output:
        txt = outdir + "/{r}/cutesv/quantify_lite/{source}.tsv"
    shell:
        """
        gzip -d -c {input.txt} \
            | awk -v OFS='\\t' '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}' > {output.txt}
        """

rule quantify_psudobulk_sv:
    input:
        vcf = outdir + "/sv/concat/cuteSV.vcf.gz",
        bam = outdir + "/{r}/merged.{hp}.bam"
    output:
        txt1 = temp(outdir + "/{r}/quantify_psudobulk_sv.{hp}.tsv"),
        txt2 = outdir + "/{r}/quantify_psudobulk_sv.{hp}.tsv.gz",
        txt3 = outdir + "/{r}/quantify_psudobulk_sv.{hp}_lite.tsv"
    threads:
        24
    shell:
        """
        ../6_SNV_SV_Comparison/scripts/quantify_sv.p.py {input.vcf} {threads} {input.bam} {output.txt1}
        pigz -p {threads} -c {output.txt1} > {output.txt2}
        awk -v OFS='\\t' '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}' {output.txt1} > {output.txt3}
        """
