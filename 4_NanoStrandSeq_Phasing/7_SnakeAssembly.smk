#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
OUTDIR = ROOT_DIR
HPS = ["hp1", "hp2"]
RS = ["round1", "round2"]

rule all:
    input:
        # expand(OUTDIR + "/{r}/bams/{chrom}.{hp}.bam", r=RS, chrom=CHROMS, hp=HPS),
        # expand(OUTDIR + "/{r}/bams/{chrom}.{hp}.flagstat", r=RS, chrom=CHROMS, hp=HPS),
        expand(OUTDIR + "/{r}/merged.bam", r=["round2"]),
        # expand(OUTDIR + "/{r}/merged.{hp}.bam", r=RS, hp=HPS),
        # expand(OUTDIR + "/{r}/quantify_psudobulk_sv.{hp}_lite.tsv", r=RS, hp=HPS),
        expand(OUTDIR + "/{r}/matrix/{chrom}.{hp}", r=RS, chrom=CHROMS, hp=HPS),
        #expand(OUTDIR + "/{r}/benchmark/{chrom}.{hp}", r=RS, chrom=CHROMS, hp=HPS),
        expand(OUTDIR + "/{r}/matrix2/{chrom}", r=RS, chrom=CHROMS),
        expand(OUTDIR + "/{r}/matrix2.stat/{chrom}.tsv", r=RS, chrom=CHROMS),
        expand(OUTDIR + "/{r}/matrix2.stat.plot/{chrom}", r=RS, chrom=CHROMS),
        expand(OUTDIR + "/{r}/matrix2.filtered/{chrom}.matrix.gz", r=RS, chrom=CHROMS),
        # expand(OUTDIR + "/{r}/snvs/{chrom}.vcf.gz", r=RS, chrom=CHROMS),
        # expand(OUTDIR + "/{r}/snvs/{chrom}.benchmark", r=RS, chrom=CHROMS),
        expand(OUTDIR + "/{r}/snvs.vcf.gz", r=RS),
        expand(OUTDIR + "/{r}/snvs_benchmark.json", r=RS),
        #expand(OUTDIR + "/{r}/snvs.benchmark", r=RS),
        # expand(OUTDIR + "/{r}/hets/{chrom}.bed", r=RS, chrom=CHROMS),
        expand(OUTDIR + "/{r}/hets.all.bed.gz", r=RS),
        #expand(OUTDIR + "/{r}/hets.benchmark", r=RS),
        #expand(OUTDIR + "/{r}/cuteSV/{chrom}.{hp}.vcf.gz", r=RS, chrom=CHROMS, hp=HPS),
        #expand(OUTDIR + "/{r}/cuteSV.{hp}.vcf.gz", r=RS, hp=HPS),
        #expand(OUTDIR + "/{r}/cutesv/quantify/{chrom}.{hp}.on_{hpA}.txt.gz", r=RS, chrom=CHROMS, hp=HPS, hpA=HPS),
        #expand(OUTDIR + "/{r}/cutesv/quantify_lite/{chrom}.{hp}.on_{hpA}.tsv", r=RS, chrom=CHROMS, hp=HPS, hpA=HPS),
        
## Round 1

rule make_round1_bam:
    input:
        txt = ROOT_DIR + "/clusters/clustered/cells.{chrom}.json"
    output:
        tmp1 = temp(OUTDIR + "/round1/bams/{chrom}.hp1.unsorted.bam"),
        tmp2 = temp(OUTDIR + "/round1/bams/{chrom}.hp2.unsorted.bam"),
        bam1 = OUTDIR + "/round1/bams/{chrom}.hp1.bam",
        bam2 = OUTDIR + "/round1/bams/{chrom}.hp2.bam"
    log:
        OUTDIR + "/round1/bams/{chrom}.log"
    threads:
        4
    shell:
        """(
        ./scripts/strandtools/make_haplotype_bam.py \
            {input.txt} {ROOT_DIR} {wildcards.chrom} \
            {output.tmp1} {output.tmp2}
        samtools sort -@ {threads} -o {output.bam1} {output.tmp1}
        samtools sort -@ {threads} -o {output.bam2} {output.tmp2}
        samtools index -@ {threads} {output.bam1}
        samtools index -@ {threads} {output.bam2} ) &> {log}
        """

## Round 2

rule split_haplotype:
    input:
        bam = ROOT_DIR + "/prepare/bams/{cell}.bam",
        bed = OUTDIR + "/round1/hets.all.bed.gz"
    output:
        out = directory(OUTDIR + "/round2/splited_haplotype/{cell}")
    log:
        OUTDIR + "/round2/splited_haplotype/{cell}.log"
    threads:
        4
    shell:
        """
        ./scripts/assembly/split_haplotype.v1.py {input.bam} {input.bed} {output.out} &> {log}
        """

rule make_round2_bam:
    input:
        expand(rules.split_haplotype.output.out, cell=CELLS),
    output:
        tmp1 = temp(OUTDIR + "/round2/bams/{chrom}.hp1.tmp.bam"),
        tmp2 = temp(OUTDIR + "/round2/bams/{chrom}.hp2.tmp.bam"),
        bam1 = OUTDIR + "/round2/bams/{chrom}.hp1.bam",
        bam2 = OUTDIR + "/round2/bams/{chrom}.hp2.bam"
    log:
        OUTDIR + "/round2/bams/{chrom}.log"
    threads:
        4
    shell:
        """(
        ./scripts/strandtools/make_haplotype_bam_r2.py \
            {CONF} {wildcards.chrom} {output.tmp1} {output.tmp2}
        samtools sort -@ {threads} -o {output.bam1} {output.tmp1}
        samtools sort -@ {threads} -o {output.bam2} {output.tmp2}
        samtools index -@ {threads} {output.bam1}
        samtools index -@ {threads} {output.bam2} ) &> {log}
        """

## Common steps

rule merge_bams:
    input:
        bams1 = [OUTDIR + "/{round}/bams/%s.hp1.bam" % c for c in CHROMS],
        bams2 = [OUTDIR + "/{round}/bams/%s.hp2.bam" % c for c in CHROMS]
    output:
        bam = OUTDIR + "/{round}/merged.bam"
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
        bams = [OUTDIR + "/{round}/bams/%s.{hp}.bam" % c for c in CHROMS]
    output:
        bam  = OUTDIR + "/{round}/merged.{hp}.bam"
    threads:
        8
    shell:
        """
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam}
        """

rule generate_base_matrix:
    input:
        bam = OUTDIR + "/{r}/bams/{chrom}.{hp}.bam",
        fa = GENOME_FASTA
    output:
        mtx = directory(OUTDIR + "/{r}/matrix/{chrom}.{hp}")
    log:
        OUTDIR + "/{r}/matrix/{chrom}.{hp}.log"
    threads:
        24
    shell:
        """
        ./scripts/assembly/generate_base_matrix.py \
            {input.bam} {input.fa} {wildcards.chrom} {threads} {output.mtx} &> {log}
        """

rule merge_base_matrix:
    input:
        fa = GENOME_FASTA,
        bed = SNP_BED,
        vcf = SNP_VCF,
        mtxdir1 = OUTDIR + "/{r}/matrix/{chrom}.hp1",
        mtxdir2 = OUTDIR + "/{r}/matrix/{chrom}.hp2"
    output:
        directory(OUTDIR + "/{r}/matrix2/{chrom}")
    log:
        OUTDIR + "/{r}/matrix2/{chrom}.log"
    threads:
        12
    shell:
        """
        ./scripts/assembly/merge_base_matrix.py {input.fa} {input.bed} {input.vcf} \
            {input.mtxdir1} {input.mtxdir2} {wildcards.chrom} {threads} {output} &> {log}
        """

rule stat_matrix2:
    input:
        OUTDIR + "/{r}/matrix2/{chrom}"
    output:
        tsv = OUTDIR + "/{r}/matrix2.stat/{chrom}.tsv"
    log:
        OUTDIR + "/{r}/matrix2.stat/{chrom}.log"
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
        directory(OUTDIR + "/{r}/matrix2.stat.plot/{chrom}")
    shell:
        """
        ./scripts/assembly/plot_stat_matrix2.py {input.tsv} {output}
        """

rule filter_matrix2:
    input:
        OUTDIR + "/{r}/matrix2/{chrom}"
    output:
        txt = OUTDIR + "/{r}/matrix2.filtered/{chrom}.matrix.gz"
    shell:
        """
        zcat {input}/*.matrix.gz \
            | awk '$6!="."&&$7!="."&&$6!="-"&&$7!="-"&&($6!=$2||$7!=$2)' \
            | sort -k1,1n | gzip -c > {output.txt}
        """

# VCF

rule make_chrom_snvs:
    input:
        mtx = OUTDIR + "/{r}/matrix2.filtered/{chrom}.matrix.gz",
        vcf = ROOT_DIR + "/snv/concat/nanocaller.vcf.gz",
        bed = ROOT_DIR + "/inversions/inversions.bed.gz",
        txt = GENOME_SIZES
    output:
        tmp = temp(OUTDIR + "/{r}/snvs/{chrom}.vcf"),
        vcf = OUTDIR + "/{r}/snvs/{chrom}.vcf.gz"
    log:
        OUTDIR + "/{r}/snvs/{chrom}.log"
    shell:
        """(
        ./scripts/assembly/call_snps.py {input.mtx} {input.vcf} {input.bed} {input.txt} {wildcards.chrom} {output.tmp}
        bgzip -c {output.tmp} > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule concat_chrom_vcfs:
    input:
        vcfs = [OUTDIR + "/{r}/snvs/%s.vcf.gz" % c for c in CHROMS]
    output:
        vcf = OUTDIR + "/{r}/snvs.vcf.gz"
    log:
        OUTDIR + "/{r}/snvs.log"
    shell:
        """(
        bcftools concat -a {input.vcfs} | bcftools sort | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule benchmark_snps:
    input:
        vcf1 = SNP_VCF,
        vcf2 = OUTDIR + "/{r}/snvs.vcf.gz",
        bed = SNP_BED
    output:
        txt = OUTDIR + "/{r}/snvs_benchmark.json"
    log:
        OUTDIR + "/{r}/snvs_benchmark.log"
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
        bed = ROOT_DIR + "/snv/concat/nanocaller.vcf.gz",
        mtxdir = OUTDIR + "/{r}/matrix2/{chrom}"
    output:
        bed = OUTDIR + "/{r}/hets/{chrom}.bed",
        bed_gz = OUTDIR + "/{r}/hets/{chrom}.bed.gz"
    log:
        OUTDIR + "/{r}/hets/{chrom}.log"
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
        beds = [OUTDIR + "/{r}/hets/%s.bed" % c for c in CHROMS]
    output:
        bed = OUTDIR + "/{r}/hets.all.bed.gz"
    shell:
        """
        cat {input.beds} | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """

# SV

rule cutesv:
    input:
        bam = OUTDIR + "/{r}/bams/{chrom}.{hp}.bam",
        fa = GENOME_FASTA
    output:
        wd = temp(directory(OUTDIR + "/{r}/cuteSV/{chrom}.{hp}.wd")),
        vcf = temp(OUTDIR + "/{r}/cuteSV/{chrom}.{hp}.vcf"),
        vcf2 = OUTDIR + "/{r}/cuteSV/{chrom}.{hp}.vcf.gz"
    log:
        OUTDIR + "/{r}/cuteSV/{chrom}.{hp}.log"
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
        vcfs = [OUTDIR + "/{r}/cuteSV/%s.{hp}.vcf.gz" % c for c in CHROMS]
    output:
        vcf = OUTDIR + "/{r}/cuteSV.{hp}.vcf.gz"
    shell:
        """
        bcftools concat -a {input.vcfs} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule quantify_sv:
    input:
        vcf = OUTDIR + "/{r}/cuteSV/{chrom}.{hp}.vcf.gz",
        bam = OUTDIR + "/{r}/bams/{chrom}.{hpA}.bam"
    output:
        tmp = temp(OUTDIR + "/{r}/cutesv/quantify/{chrom}.{hp}.on_{hpA}.txt"),
        txt = OUTDIR + "/{r}/cutesv/quantify/{chrom}.{hp}.on_{hpA}.txt.gz"
    threads:
        THREADS
    shell:
        """
        ../3_NanoStrandSeq_PseudoBulk/scripts/quantify_sv.p.py {input.vcf} {threads} {input.bam} {output.tmp}
        awk 'NR==1||$1=="{wildcards.chrom}"' {output.tmp} | pigz -p {threads} -c > {output.txt}
        """

rule lite_quantify:
    input:
        txt = OUTDIR + "/{r}/cutesv/quantify/{source}.txt.gz"
    output:
        txt = OUTDIR + "/{r}/cutesv/quantify_lite/{source}.tsv"
    shell:
        """
        gzip -d -c {input.txt} \
            | awk -v OFS='\\t' '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}' > {output.txt}
        """

rule quantify_psudobulk_sv:
    input:
        vcf = OUTDIR + "/sv/concat/cuteSV.vcf.gz",
        bam = OUTDIR + "/{r}/merged.{hp}.bam"
    output:
        txt1 = temp(OUTDIR + "/{r}/quantify_psudobulk_sv.{hp}.tsv"),
        txt2 = OUTDIR + "/{r}/quantify_psudobulk_sv.{hp}.tsv.gz",
        txt3 = OUTDIR + "/{r}/quantify_psudobulk_sv.{hp}_lite.tsv"
    threads:
        24
    shell:
        """
        ../3_NanoStrandSeq_PseudoBulk/scripts/quantify_sv.p.py {input.vcf} {threads} {input.bam} {output.txt1}
        pigz -p {threads} -c {output.txt1} > {output.txt2}
        awk -v OFS='\\t' '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}' {output.txt1} > {output.txt3}
        """
