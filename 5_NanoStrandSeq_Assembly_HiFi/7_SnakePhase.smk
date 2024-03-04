#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = "results/phase"
# chroms = ["cluster_26"]
chroms = clusters
hps = ["hp1", "hp2"]

rule all:
    input:
        expand(outdir + "/wc/{cell}", cell=nano_strand_seq_cells),
        outdir + "/anchors.bed.gz",
        expand(outdir + "/assigned/{cell}.bed.gz", cell=nano_strand_seq_cells),
        outdir + "/clustered",
        expand(outdir + "/round1/bams/{chrom}.{hp}.bam", chrom=clusters, hp=hps),
        expand(outdir + "/round1/matrix/{chrom}.{hp}", chrom=clusters, hp=hps),
        expand(outdir + "/round1/matrix2/{chrom}", chrom=chroms),
        expand(outdir + "/round1/matrix2.stat/{chrom}.tsv", chrom=chroms),
        expand(outdir + "/round1/matrix2.stat.plot/{chrom}", chrom=chroms),
        expand(outdir + "/round1/matrix2.filtered/{chrom}.matrix.gz", chrom=chroms),
        expand(outdir + "/round1/snvs/{chrom}.vcf.gz", chrom=chroms),
        expand(outdir + "/whatshap/unphase/{chrom}.vcf.gz", chrom=chroms),
        expand(outdir + "/whatshap/phased_ccs/{chrom}.vcf.gz", chrom=chroms),
        outdir + "/whatshap/phased_ccs.vcf.gz",
        # expand(outdir + "/whatshap/ccs_haplotag/{chrom}.bam", chrom=chroms),


rule fetch_wc_region:
    input:
        bam = "results/ss2cluster/tgs/mark_duplicate/{cell}.bam"
    output:
        out = directory(outdir + "/wc/{cell}")
    log:
        outdir + "/wc/{cell}.log"
    shell:
        """
        sstools FetchWCRegion {input.bam} {output.out} &> {log}
        """

rule get_anchors:
    input:
        vcf = "results/ccs2cluster/nanocaller.vcf.gz"
    output:
        bed = outdir + "/anchors.bed.gz"
    shell:
        """
        zcat {input.vcf} | grep -v '#' | grep -v '1/1' | awk 'length($4)==1 && length($5)==1' \
            | awk -v FS='\\t' -v OFS='\\t' '{{print $1,$2-1,$2,$4"/"$5}}' \
            | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """

rule assign_anchor:
    input:
        bam = "results/ss2cluster/tgs/mark_duplicate/{cell}.bam",
        bed = rules.get_anchors.output.bed,
        txt = outdir + "/wc/{cell}/duplicate_set_names.json"
    output:
        bed = outdir + "/assigned/{cell}.bed.gz"
    log:
        outdir + "/assigned/{cell}.log"
    shell:
        """(
        ../1_NanoStrandseq/scripts/strandtools/assign_anchor.py {input.bam} {input.bed} {input.txt} \
            | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
        tabix -p bed {output.bed} ) &> {log}
        """

rule clustering:
    input:
        beds = expand(outdir + "/assigned/{cell}.bed.gz", cell=nano_strand_seq_cells),
    output:
        txt = outdir + "/clustered.config",
        out = directory(outdir + "/clustered")
    log:
        outdir + "/clustered.log"
    params:
        chroms = ",".join(chroms)
    shell:
        """
        for f in {input.beds}; do
            cell=`basename $f .bed.gz`
            echo -e "${{cell}}\\t${{f}}"
        done > {output.txt}
        ../1_NanoStrandseq/scripts/strandtools/cluster_cells.v2.py {output.txt} {params.chroms} {output.out} &> {log}
        """


rule make_round1_bam:
    input:
        # txt = assembly_dir + "/clusters/clustered/cells.{chrom}.json"
        txtdir = outdir + "/clustered"
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
        ../1_NanoStrandseq/scripts/strandtools/make_haplotype_bam.v2.py \
            {input.txtdir}/cells.{wildcards.chrom}.json None {wildcards.chrom} \
            {output.tmp1} {output.tmp2}
        samtools sort -@ {threads} -o {output.bam1} {output.tmp1}
        samtools sort -@ {threads} -o {output.bam2} {output.tmp2}
        samtools index -@ {threads} {output.bam1}
        samtools index -@ {threads} {output.bam2} ) &> {log}
        """

rule generate_base_matrix:
    input:
        bam = outdir + "/{r}/bams/{chrom}.{hp}.bam",
        fa = cluster_fasta
    output:
        mtx = directory(outdir + "/{r}/matrix/{chrom}.{hp}")
    log:
        outdir + "/{r}/matrix/{chrom}.{hp}.log"
    threads:
        24
    shell:
        """
        ../1_NanoStrandseq/scripts/assembly/generate_base_matrix.py \
            {input.bam} {input.fa} {wildcards.chrom} {threads} {output.mtx} &> {log}
        """

rule merge_base_matrix:
    input:
        fa = cluster_fasta,
        # bed = lambda wildcards: BENCHMARKS[cellline]["BED"],
        # vcf = lambda wildcards: BENCHMARKS[cellline]["VCF"],
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
        ../1_NanoStrandseq/scripts/assembly/merge_base_matrix.v2.py {input.fa} BED VCF \
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
        ../1_NanoStrandseq/scripts/assembly/stat_matrix2.py {input} {threads} {output.tsv} &> {log}
        """

rule plot_stat_matrix2:
    input:
        tsv = rules.stat_matrix2.output.tsv
    output:
        directory(outdir + "/{r}/matrix2.stat.plot/{chrom}")
    shell:
        """
        ../1_NanoStrandseq/scripts/assembly/plot_stat_matrix2.py {input.tsv} {output}
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

rule make_chrom_snvs:
    input:
        mtx = outdir + "/{r}/matrix2.filtered/{chrom}.matrix.gz",
        vcf = "results/ccs2cluster/nanocaller/{chrom}.vcf.gz",
        # bed = assembly_dir + "/inversions/inversions.bed.gz",
        txt = cluster_sizes
    output:
        tmp = temp(outdir + "/{r}/snvs/{chrom}.vcf"),
        vcf = outdir + "/{r}/snvs/{chrom}.vcf.gz"
    log:
        outdir + "/{r}/snvs/{chrom}.log"
    shell:
        """(
        ../1_NanoStrandseq/scripts/assembly/call_snps.v2.py {input.mtx} {input.vcf} BED {input.txt} {wildcards.chrom} {output.tmp}
        bgzip -c {output.tmp} > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule unphase:
    input:
        vcf = "results/ccs2cluster/nanocaller/{chrom}.vcf.gz"
    output:
        vcf = temp(outdir + "/whatshap/unphase/{chrom}.vcf"),
        vcf_gz = outdir + "/whatshap/unphase/{chrom}.vcf.gz"
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

rule phase_nss:
    input:
        fa = cluster_fasta,
        vcf1 = outdir + "/whatshap/unphase/{chrom}.vcf.gz",
        vcf2 = outdir + "/round1/snvs/{chrom}.vcf.gz",
        bam = "results/ccs2cluster/minimap2/mapped_merged_filtered.bam"
    output:
        vcf = outdir + "/whatshap/phased_ccs/{chrom}.vcf.gz"
    log:
        outdir + "/whatshap/phased_ccs/{chrom}.log"
    shell:
        """(
        whatshap phase --ignore-read-groups -r {input.fa} --chromosome {wildcards.chrom} \
            {input.vcf1} {input.vcf2} {input.bam} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf} ) &> {log}
        """

rule merge_vcf:
    input:
        vcfs = expand(outdir + "/whatshap/phased_ccs/{chrom}.vcf.gz", chrom=chroms)
    output:
        vcf = outdir + "/whatshap/phased_ccs.vcf.gz"
    shell:
        """
        bcftools concat -a {input.vcfs} | bcftools sort | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule haplotag:
    input:
        fa = cluster_fasta,
        bam = "results/ccs2cluster/minimap2/mapped_merged_filtered.bam",
        vcf = outdir + "/whatshap/phased_ccs/{chrom}.vcf.gz"
    output:
        bam = outdir + "/whatshap/ccs_haplotag/{chrom}.bam",
        txt = outdir + "/whatshap/ccs_haplotag/{chrom}.txt"
    log:
        outdir + "/whatshap/ccs_haplotag/{chrom}.log"
    shell:
        """
        whatshap haplotag \
            --ignore-read-groups \
            --regions {wildcards.chrom} \
            --output {output.bam} \
            --reference {input.fa} \
            --output-haplotag-list {output.txt} \
            {input.vcf} {input.bam} &> {log}
        """