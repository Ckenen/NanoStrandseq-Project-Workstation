#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
OUTDIR = ROOT_DIR + "/prepare"
# chroms = ["chr22"]

rule all:
    input:
        OUTDIR + "/config.json",
        expand(OUTDIR + "/bams/{cell}.bam", cell=CELLS),
        expand(OUTDIR + "/stat_bin_reads/{cell}.RmDup0.tsv", cell=CELLS),
        expand(OUTDIR + "/plot_bin_reads/{cell}.RmDup0.pdf", cell=CELLS),
        # OUTDIR + "/all_cells.barplot.RmDup0.pdf",
        OUTDIR + "/all_cells.all_chroms.bam",
        OUTDIR + "/all_cells.all_chroms.flagstat",
        OUTDIR + "/all_cells.all_chroms.bw",
        OUTDIR + "/all_cells.all_chroms.rmdup.bam",
        OUTDIR + "/all_cells.all_chroms.rmdup.flagstat",
        OUTDIR + "/all_cells.all_chroms.rmdup.bw",
        OUTDIR + "/all_cells.all_chroms.raw.bam",
        OUTDIR + "/all_cells.all_chroms.raw.flagstat",
        OUTDIR + "/all_cells.all_chroms.raw.bw",
        # OUTDIR + "/breakpoints.bedGraph",
        # OUTDIR + "/all_cells.all_chroms.raw.cutesv.vcf.gz",

rule copy_config:
    input:
        txt = CONF
    output:
        txt = OUTDIR + "/config.json"
    shell:
        """
        cp {input.txt} {output.txt}
        """

rule link_bams:
    input:
        bam = lambda wildcards: get_bam(wildcards.cell)
    output:
        bam = OUTDIR + "/bams/{cell}.bam"
    shell:
        """
        bam=`readlink -f {input.bam}`
        ln -s -f ${{bam}} {output.bam}
        ln -s -f ${{bam}}.bai {output.bam}.bai
        """

rule stat_bin_reads:
    input:
        bam = OUTDIR + "/bams/{cell}.bam"
    output:
        tsv1 = OUTDIR + "/stat_bin_reads/{cell}.RmDup0.tsv",
        tsv2 = OUTDIR + "/stat_bin_reads/{cell}.RmDup1.tsv"
    threads:
        8
    shell:
        """
        sstools StatBinRead -w 1000000 -t {threads} {input.bam} {output.tsv1}
        sstools StatBinRead -w 1000000 -t {threads} -d {input.bam} {output.tsv2}
        """

rule plot_bin_reads:
    input:
        tsv1 = rules.stat_bin_reads.output.tsv1,
        tsv2 = rules.stat_bin_reads.output.tsv2
    output:
        pdf1 = OUTDIR + "/plot_bin_reads/{cell}.RmDup0.pdf",
        pdf2 = OUTDIR + "/plot_bin_reads/{cell}.RmDup1.pdf"
    shell:
        """
        sstools PlotBinRead -t {input.tsv1} {output.pdf1}
        sstools PlotBinRead -t {input.tsv2} {output.pdf2}
        """

rule merge_pdfs:
    input:
        pdfs1 = expand(OUTDIR + "/plot_bin_reads/{cell}.RmDup0.pdf", cell=CELLS),
        pdfs2 = expand(OUTDIR + "/plot_bin_reads/{cell}.RmDup1.pdf", cell=CELLS)
    output:
        pdf1 = OUTDIR + "/all_cells.barplot.RmDup0.pdf",
        pdf2 = OUTDIR + "/all_cells.barplot.RmDup1.pdf"
    shell:
        """
        merge_pdf.py {input.pdfs1} {output.pdf1}
        merge_pdf.py {input.pdfs2} {output.pdf2}
        """

rule merge_bams:
    input:
        bams = expand(OUTDIR + "/bams/{cell}.bam", cell=CELLS)
    output:
        bam = OUTDIR + "/all_cells.all_chroms.bam"
    threads:
        THREADS
    shell:
        """
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam}
        """

rule make_rmdup_bam:
    input:
        bam = OUTDIR + "/all_cells.all_chroms.bam"
    output:
        bam = OUTDIR + "/all_cells.all_chroms.rmdup.bam"
    threads:
        THREADS
    shell:
        """
        samtools view -@ {threads} -F 1024 -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam}
        """

rule merge_raw_bams:
    input:
        bams = [get_raw_bam(cell) for cell in CELLS]
    output:
        bam = OUTDIR + "/all_cells.all_chroms.raw.bam"
    threads:
        THREADS
    shell:
        """
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam}
        """

rule parse_breakpoints:
    input:
        bam = rules.merge_bams.output.bam
    output:
        bg = OUTDIR + "/breakpoints.bedGraph"
    log:
        OUTDIR + "/breakpoints.log"
    threads:
        48
    shell:
        """
        ./scripts/parse_breakpoints.py {input.bam} {threads} {output.bg} &> {log}
        """

# Common rules

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    shell:
        """
        samtools flagstat {input} > {output}
        """

rule bam_to_bw:
    input:
        bam = "{prefix}.bam",
        txt = GENOME_SIZES
    output:
        td = temp(directory("{prefix}.BAM_TO_BW_SORT_TMP")),
        bed = temp("{prefix}.bed"),
        bg = temp("{prefix}.bedGraph"),
        bed2 = "{prefix}.bed.gz",
        bw = "{prefix}.bw"
    threads:
        8
    shell:
        """
        mkdir -p {output.td}
        bedtools bamtobed -i {input.bam} \
            | awk -v FS='\\t' -v OFS='\\t' '{{if($6=="+"){{c="107,137,138"}}else{{c="248,173,97"}};print $1,$2,$3,$4,$5,$6,$2,$3,c}}' \
            | sort --parallel {threads} -T {output.td} -k1,1 -k2,2n -k3,3n > {output.bed}
        bgzip -c {output.bed} > {output.bed2}
        tabix -p bed {output.bed2}
        bedtools genomecov -bg -i {output.bed} -g {input.txt} | sort --parallel {threads} -T {output.td} -k1,1 -k2,2n > {output.bg}
        bedGraphToBigWig {output.bg} {input.txt} {output.bw}
        """

# rule cutesv:
#     input:
#         fasta = lambda wildcards: GENOMES[species]["GENOME_FASTA"],
#         bam = OUTDIR + "/all_cells.all_chroms.raw.bam"
#     output:
#         wd = temp(directory(OUTDIR + "/all_cells.all_chroms.raw.cutesv.wd")),
#         vcf = temp(OUTDIR + "/all_cells.all_chroms.raw.cutesv.vcf"),
#         vcf2 = OUTDIR + "/all_cells.all_chroms.raw.cutesv.vcf.gz"
#     log:
#         OUTDIR + "/all_cells.all_chroms.raw.cutesv.log"
#     threads:
#         48 # threads
#     shell:
#         """(
#         mkdir -p {output.wd}
#         cuteSV --max_cluster_bias_INS 100 \
#             --diff_ratio_merging_INS 0.3 \
#             --max_cluster_bias_DEL 100 \
#             --diff_ratio_merging_DEL 0.3 \
#             --threads {threads} \
#             --sample cuteSV \
#             --retain_work_dir \
#             --report_readid \
#             --min_support 1 \
#             --min_size 50 \
#             --genotype \
#             {input.bam} {input.fasta} {output.vcf} {output.wd} 
#         bgzip -c {output.vcf} > {output.vcf2}
#         tabix -p vcf {output.vcf2} ) &> {log}
#         """