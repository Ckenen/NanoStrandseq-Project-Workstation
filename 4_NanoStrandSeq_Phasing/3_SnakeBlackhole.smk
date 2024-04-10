#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
OUTDIR = ROOT_DIR + "/blackhole"

rule all:
    input:
        OUTDIR + "/bin_reads.bed.gz",
        OUTDIR + "/blacklist.bed.gz",
        OUTDIR + "/whitelist.bed.gz"

rule cal_bin_reads:
    input:
        bam = ROOT_DIR + "/prepare/all_cells.all_chroms.bam"
    output:
        bed = OUTDIR + "/bin_reads.bed",
        bed2 = OUTDIR + "/bin_reads.bed.gz"
    threads:
        24
    shell:
        """
        ./scripts/strandtools/cal_bin_reads.py {input.bam} {threads} {output.bed}
        sort -k1,1 -k2,2n {output.bed} | bgzip -c > {output.bed2}
        tabix -p bed {output.bed2}
        """

rule call_blacklist:
    input:
        bed = rules.cal_bin_reads.output.bed
    output:
        bed = OUTDIR + "/blacklist.bed",
        bed2 = OUTDIR + "/blacklist.bed.gz"
    shell:
        """
        ./scripts/assembly/call_blacklist.py {input.bed} {output.bed}
        bedtools merge -i {output.bed} | sort -k1,1 -k2,2n | bgzip -c > {output.bed2}
        tabix -p bed {output.bed2}
        """

rule get_whitelist:
    input:
        bed = rules.call_blacklist.output.bed2,
        txt = GENOME_SIZES
    output:
        txt = temp(OUTDIR + "/whitelist.sizes"),
        bed = OUTDIR + "/whitelist.bed.gz"
    shell:
        """
        sort -k1,1 {input.txt} > {output.txt}
        bedtools complement -i {input.bed} -g {output.txt} | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """