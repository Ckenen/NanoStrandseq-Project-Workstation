#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "data/datasets"
outdir = "results/demux"

rule all:
    input:
        expand(outdir + "/barcodes/{run}.1.fa", run=runs),
        expand(outdir + "/fbilr/{run}.matrix.gz", run=runs),
        expand(outdir + "/fbilr/{run}.stats.tsv.gz", run=runs),
        # expand(outdir + "/splitted/{run}", run=runs),
        # expand(outdir + "/combined/{run_cell}.fastq.gz", run_cell=run_cells),
        expand(outdir + "/trimmed/{run_cell}", run_cell=run_cells),

rule get_barcodes:
    input:
        fa = config["barcodes"]
    output:
        fa1 = outdir + "/barcodes/{run}.1.fa",
        fa2 = outdir + "/barcodes/{run}.2.fa",
        tsv = outdir + "/barcodes/{run}.tsv"
    run:
        import subprocess
        d = dat[dat["Run"] == wildcards.run]
        bcs1 = ["Bar%d" % bc for bc in sorted(set(d["Barcode.1st"]))]
        bcs2 = ["Bar%d" % bc for bc in sorted(set(d["Barcode.2nd"]))]
        cmd1 = "samtools faidx %s %s > %s" % (input.fa, " ".join(bcs1), output.fa1)
        cmd2 = "samtools faidx %s %s > %s" % (input.fa, " ".join(bcs2), output.fa2)
        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
        with open(output.tsv, "w+") as fw:
            for cell, bc1, bc2 in d[["Cell", "Barcode.1st", "Barcode.2nd"]].values:
                fw.write("%s\tBar%s\tBar%s\n" % (cell, bc1, bc2))

rule fbilr:
    input:
        fq = indir + "/{run}.fastq.gz",
        fa1 = rules.get_barcodes.output.fa1,
        fa2 = rules.get_barcodes.output.fa2
    output:
        mtx = outdir + "/fbilr/{run}.matrix.gz"
    log:
        outdir + "/fbilr/{run}.log"
    threads:
        12
    shell:
        """
        fbilr -t {threads} -w 200 -b {input.fa1},{input.fa2} {input.fq} 2> {log} \
            | pigz -p {threads} -c > {output.mtx}
        """

rule report_matrix_summary:
    input:
        mtx = rules.fbilr.output.mtx
    output:
        tsv = outdir + "/fbilr/{run}.stats.tsv.gz"
    shell:
        """
        zcat {input.mtx} | awk '$2>=400' | awk -v OFS=',' '{{print $3,$4,$5,$8,$9,$10,$11,$14}}' \
            | sort | uniq -c | awk -v OFS=',' '{{print $2,$1}}' \
            | sed 's/,/\\t/g' | gzip -c > {output.tsv}
        """

rule split_reads:
    input:
        fq = indir + "/{run}.fastq.gz",
        mtx = rules.fbilr.output.mtx,
        tsv = rules.get_barcodes.output.tsv
    output:
        out = directory(outdir + "/splitted/{run}")
    log:
        outdir + "/splitted/{run}.log"
    threads:
        4
    shell:
        """
        nss_split_reads.py -e 6 -l 400 {input.fq} {input.mtx} {input.tsv} {output.out} &> {log}
        """

rule combine_reads:
    input:
        fqs = rules.split_reads.output.out
    output:
        fq = outdir + "/combined/{run}/{cell}.fastq.gz"
    threads:
        4
    shell:
        """
        fq1="{input.fqs}/fastqs/{wildcards.cell}_F.fastq"
        fq2="{input.fqs}/fastqs/{wildcards.cell}_R.fastq"
        ( cat $fq1; cat $fq2 | reverse_fastq.py ) | pigz -p {threads} -c > {output.fq}
        """

# Trimmed

rule trim_reads:
    input:
        fq = outdir + "/combined/{run}/{cell}.fastq.gz"
    output:
        out = directory(outdir + "/trimmed/{run}/{cell}")
    log:
        outdir + "/trimmed/{run}/{cell}.log"
    shell:
        """
        nss_trim_reads.py {input.fq} {output.out} &> {log}
        """
