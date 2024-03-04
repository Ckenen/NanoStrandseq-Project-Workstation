#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
outdir = assembly_dir + "/inversions.v2"

rule all:
    input:
        outdir + "/regions/filelist.tsv",
        outdir + "/regions/all_cc_regions.bed.gz",
        outdir + "/composites/bamlist.tsv",
        outdir + "/composites/all_chroms.bed.gz",
        outdir + "/composites/all_chroms.+.bw",
        outdir + "/composites/all_chroms.-.bw",
        outdir + "/inversions.bed.gz"

rule make_count_list:
    input:
        tsv_list = [assembly_dir + "/prepare/stat_bin_reads/%s.RmDup1.tsv" % cell for cell in cells]
    output:
        tsv = outdir + "/regions/filelist.tsv"
    run:
        with open(output.tsv, "w+") as fw:
            for path in input.tsv_list:
                cell = os.path.basename(path)[:-11]
                fw.write("%s\t%s\n" % (cell, path))
    
rule fetch_cc_regions:
    input:
        tsv = rules.make_count_list.output.tsv
    output:
        out = directory(outdir + "/regions/all_cc_regions.outdir"),
        bed = outdir + "/regions/all_cc_regions.bed.gz"
    log:
        outdir + "/regions/all_cc_regions.log"
    threads:
        threads
    shell:
        """(
        sstools FetchCCRegion -t {threads} {input.tsv} {output.out}
        sort -k1,1 -k2,2n {output.out}/all_cc_regions.bed | bgzip -c > {output.bed}
        tabix -p bed {output.bed} ) &> {log}
        """

rule make_bam_list:
    input:
        bam_list = [assembly_dir + "/prepare/bams/%s.bam" % cell for cell in cells]
    output:
        tsv = outdir + "/composites/bamlist.tsv"
    run:
        with open(output.tsv, "w+") as fw:
            for path in input.bam_list:
                cell = os.path.basename(path)[:-4]
                fw.write("%s\t%s\n" % (cell, path))

rule make_cc_composite:
    input:
        bed = rules.fetch_cc_regions.output.bed,
        tsv = rules.make_bam_list.output.tsv
    output:
        out = temp(directory(outdir + "/composites/chroms")),
        bed = outdir + "/composites/all_chroms.bed.gz"
    log:
        outdir + "/composites/all_chroms.log"
    threads:
        8
    shell:
        """(
        sstools MakeCCComposite -t {threads} {input.bed} {input.tsv} {output.out}
        cat {output.out}/*.bed | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
        tabix -p bed {output.bed} ) &> {log}
        """

rule make_bigwig:
    input:
        bed = "{prefix}.bed.gz",
        txt = lambda wildcards: GENOMES[species]["GENOME_SIZE"]
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
        bed = rules.make_cc_composite.output.bed
    output:
        bed1 = temp(outdir + "/inversions.bed"),
        bed2 = outdir + "/inversions.bed.gz"
    log:
        outdir + "/inversions.log"
    threads:
        24
    shell:
        """(
        sstools CallInversion -t {threads} {input.bed} {output.bed1}
        sort -k1,1 -k2,2n {output.bed1} | bgzip -c > {output.bed2}
        tabix -p bed {output.bed2} ) &> {log}
        """