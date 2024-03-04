#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/demux/trimmed"
outdir = "results/mapping"
# run_cells = run_cells[:1]

rule all:
    input:
        expand(outdir + "/minimap2/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/minimap2/{run_cell}.flagstat", run_cell=run_cells),
        # expand(outdir + "/filtered/{run_cell}.bam", run_cell=run_cells),
        # expand(outdir + "/mark_region/{run_cell}.bam", run_cell=run_cells),
        # expand(outdir + "/mark_parental/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/mark_duplicate/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/mark_duplicate/{run_cell}.flagstat", run_cell=run_cells),
        expand(outdir + "/remove_duplicate/{run_cell}.bam", run_cell=run_cells),
        expand(outdir + "/remove_duplicate/{run_cell}.flagstat", run_cell=run_cells),
        outdir + "/final_reads.%s.tsv" % config["name"]

rule minimap2:
    input:
        fq = indir + "/{run}/{cell}/trimmed.fastq.gz",
        mmi = config["mmi"]
    output:
        bam = outdir + "/minimap2/{run}/{cell}.bam"
    log:
        outdir + "/minimap2/{run}/{cell}.log"
    params:
        rg = "@RG\\tID:{cell}\\tLB:{cell}\\tSM:{cell}"
    threads:
        threads
    shell:
        """(
        minimap2 -ax map-ont --MD -t {threads} -R '{params.rg}' --secondary=no {input.mmi} {input.fq} \
            | samtools view -@ {threads} -u -F 4 - \
            | samtools sort -@ {threads} -T {output.bam}_TMP -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule filter_bam:
    input:
        bam = rules.minimap2.output.bam
    output:
        bam = outdir + "/filtered/{run}/{cell}.bam"
    log:
        outdir + "/filtered/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools FilterBam -n '^chr([0-9]+|[XY])$' -q 30 {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_region:
    input:
        bam = rules.filter_bam.output.bam,
        bed = config["benchmark_bed"]
    output:
        bam = outdir + "/mark_region/{run}/{cell}.bam"
    log:
        outdir + "/mark_region/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools MarkRegion -n XH {input.bam} {input.bed} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_parental:
    input:
        bam = rules.mark_region.output.bam,
        vcf = config["benchmark_vcf"]
    output:
        bam = outdir + "/mark_parental/{run}/{cell}.bam"
    log:
        outdir + "/mark_parental/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools MarkHaplotype -n XP {input.bam} {input.vcf} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """

rule mark_duplicate:
    input:
        bam = rules.mark_parental.output.bam
    output:
        bam = outdir + "/mark_duplicate/{run}/{cell}.bam"
    log:
        outdir + "/mark_duplicate/{run}/{cell}.log"
    threads:
        4
    shell:
        """
        sstools MarkDuplicate -d 20 -s {input.bam} {output.bam} &> {log}
        samtools index -@ {threads} {output.bam}
        """


rule remove_duplicate:
    input:
        bam = rules.mark_duplicate.output.bam
    output:
        bam = outdir + "/remove_duplicate/{run}/{cell}.bam"
    threads:
        4
    shell:
        """
        samtools view -@ {threads} -F 3072 -b {input.bam} > {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule report_final_reads:
    input:
        txts = expand(outdir + "/mark_duplicate/{run_cell}.flagstat", run_cell=run_cells)
    output:
        tsv = outdir + "/final_reads.%s.tsv" % config["name"]
    run:
        with open(output.tsv, "w+") as fw:
            fw.write("Cell\tReads\tDupReads\tUniqReads\n")
            for path in sorted(input.txts):
                cell = path.split("/")[-1][:-9]
                total, uniq, dup = 0, 0, 0
                for line in open(path):
                    if "primary mapped" in line:
                        total = int(line.split()[0])
                    if "primary duplicates" in line:
                        dup = int(line.split()[0])
                line = "\t".join(map(str, [cell, total, dup, total - dup]))
                fw.write(line + "\n")


# Common rules

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
