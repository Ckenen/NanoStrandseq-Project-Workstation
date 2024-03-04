#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = assembly_dir + "/inversions"


rule all:
    input:
        #outdir + "/regions/cc_regions.bed.gz",
        #outdir + "/composites/all_chroms.bed.gz",
        #outdir + "/composites/all_chroms.+.bw",
        #outdir + "/composites/all_chroms.-.bw",
        outdir + "/inversions.bed.gz"

rule fetch_cc_regions:
    input:
        txt = assembly_dir + "/prepare/config.json"
    output:
        out = directory(outdir + "/regions/cc_regions.outdir"),
        bed = outdir + "/regions/cc_regions.bed.gz",
        tbi = outdir + "/regions/cc_regions.bed.gz.tbi"
    log:
        outdir + "/regions/cc_regions.log"
    threads:
        threads
    shell:
        """(
        ./scripts/strandtools/fetch_cc_regions.py -p {threads} {input.txt} {output.out}
        sort -k1,1 -k2,2n {output.out}/cc_regions.bed | bgzip -c > {output.bed}
        tabix -p bed {output.bed} ) &> {log}
        """
        # """(
        # sstools FetchCCRegion -p {threads} {input.txt} {output.out}
        # sort -k1,1 -k2,2n {output.out}/cc_regions.bed | bgzip -c > {output.bed}
        # tabix -p bed {output.bed} ) &> {log}
        # """

rule make_composite:
    input:
        bed = rules.fetch_cc_regions.output.bed
    output:
        out = temp(directory(outdir + "/composites/chroms")),
        bed = outdir + "/composites/all_chroms.bed.gz",
        tbi = outdir + "/composites/all_chroms.bed.gz.tbi"
    log:
        outdir + "/composites/all_chroms.log"
    threads:
        4
    shell:
        """(
        sstools MakeCCComposite -p {threads} {input.bed} {output.out}
        cat {output.out}/*.bed | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
        tabix -p bed {output.bed} ) &> {log}
        """

rule make_bigwig:
    input:
        bed = "{prefix}.bed.gz",
        txt = lambda wildcards: get_genome_size(cells[0])
    output:
        bg = temp("{prefix}.{s}.bg"),
        bw = "{prefix}.{s}.bw"
    shell:
        """
        zcat {input.bed} \
            | bedtools genomecov -bg -strand {wildcards.s} -i - -g {input.txt} \
            | sort -k1,1 -k2,2n > {output.bg}
        bedGraphToBigWig {output.bg} {input.txt} {output.bw}
        """

rule call_inversion:
    input:
        bed = rules.make_composite.output.bed
    output:
        bed1 = temp(outdir + "/inversions.bed"),
        bed2 = outdir + "/inversions.bed.gz",
        tbi2 = outdir + "/inversions.bed.gz.tbi"
    log:
        outdir + "/inversions.log"
    threads:
        24
    shell:
        """(
        sstools CallInversion -p {threads} {input.bed} {output.bed1}
        sort -k1,1 -k2,2n {output.bed1} | bgzip -c > {output.bed2}
        tabix -p bed {output.bed2} ) &> {log}
        """