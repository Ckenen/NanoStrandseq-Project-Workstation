#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/prepare/cutadapt"
outdir = "results/mapping"

rule all:
    input:
        expand(outdir + "/bwa/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/bwa/{run_cell}.flagstat", run_cell=run_cells),
        expand(outdir + "/filtered/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/filtered/{run_cell}.flagstat", run_cell=run_cells),
        expand(outdir + "/mark_duplicates/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/mark_duplicates/{run_cell}.flagstat", run_cell=run_cells),
        expand(outdir + "/remove_duplicates/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/mark_confident_region/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/mark_parental/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/mark_parental/{run_cell}.flagstat", run_cell=run_cells),

rule bwa_mem:
    input:
        fq1 = indir + "/{run}/{cell}_1.fastq.gz",
        fq2 = indir + "/{run}/{cell}_2.fastq.gz",
        idx = lambda wildcards: GENOMES[get_species(wildcards.cell)]["GENOME_BWA"]
    output:
        bam = outdir + "/bwa/{run}/{cell}.bam"
    log:
        outdir + "/bwa/{run}/{cell}.log"
    params:
        rg = '@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}'
    threads:
        8
    shell:
        """(
        bwa mem -t {threads} -R \'{params.rg}\' {input.idx}/ref {input.fq1} {input.fq2} \
            | samtools view -F 4 -u - \
            | samtools sort -@ {threads} -T {output.bam}_TMP - > {output.bam} 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.bwa_mem.output.bam
    output:
        bam = outdir + "/filtered/{run}/{cell}.bam"
    log:
        outdir + "/filtered/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        ./scripts/filter_paired_end_bam.py {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_duplicates:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = outdir + "/mark_duplicates/{run}/{cell}.bam",
        txt = outdir + "/mark_duplicates/{run}/{cell}.metrics"
    log:
        outdir + "/mark_duplicates/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        picard MarkDuplicates --REMOVE_DUPLICATES false -I {input.bam} -M {output.txt} -O {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """


rule remove_duplicates:
    input:
        bam = rules.filter_bam.output.bam
    output:
        bam = outdir + "/remove_duplicates/{run}/{cell}.bam",
        txt = outdir + "/remove_duplicates/{run}/{cell}.metrics"
    log:
        outdir + "/remove_duplicates/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        picard MarkDuplicates --REMOVE_DUPLICATES true -I {input.bam} -M {output.txt} -O {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_confident_region:
    input:
        bam = rules.remove_duplicates.output.bam,
        bed = lambda wildcards: BENCHMARKS[get_cellline(wildcards.cell)]["BED"]
    output:
        bam = outdir + "/mark_confident_region/{run}/{cell}.bam"
    log:
        outdir + "/mark_confident_region/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools MarkRegion -f {input.bed} -n XH {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_parental:
    input:
        bam = rules.mark_confident_region.output.bam,
        vcf = lambda wildcards: BENCHMARKS[get_cellline(wildcards.cell)]["VCF"]
    output:
        bam = outdir + "/mark_parental/{run}/{cell}.bam"
    log:
        outdir + "/mark_parental/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools MarkHaplotype -n XP -p {input.vcf} {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

# rule final_reads:
#     input:
#         expand(outdir + "/mark_parental/{run_cell}.flagstat", run_cell=run_cells),
#     output:
#         tsv = outdir + "/final_reads.tsv"
#     shell:
#         """
#         echo -e "Cell\\tReads\\tDupReads\\tUniqReads" > {output.tsv}
#         for path in {input}; do
#             cell=`basename $path .flagstat`
#             reads=`cat $path | awk '{{print $1}}' | head -n 1`
#             dups=`cat $path | awk '{{print $1}}' | head -n 5 | tail -n 1`
#             echo $cell $reads $dups | awk -v OFS='\\t' '{{print $1,$2,$3,$2-$3}}'
#         done >> {output.tsv}
#         """

# Common rule

rule bam_flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """