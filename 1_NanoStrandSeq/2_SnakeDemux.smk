#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
INDIR = "data/datasets"
OUTDIR = "results/demux"

rule all:
    input:
        expand(OUTDIR + "/barcodes/{run}.1.fa", run=RUNS),
        expand(OUTDIR + "/fbilr/{run}.tsv.gz", run=RUNS),
        expand(OUTDIR + "/fbilr_stats/{run}.tsv", run=RUNS),
        expand(OUTDIR + "/splitted/{run}", run=RUNS),
        expand(OUTDIR + "/combined/{run_cell}.fastq.gz", run_cell=RUN_CELLS),
        expand(OUTDIR + "/trimmed/{run_cell}", run_cell=RUN_CELLS),

rule get_barcodes:
    input:
        fa = config["BARCODES"]
    output:
        fa1 = OUTDIR + "/barcodes/{run}.1.fa",
        fa2 = OUTDIR + "/barcodes/{run}.2.fa",
        tsv = OUTDIR + "/barcodes/{run}.tsv"
    run:
        import subprocess
        d = DAT[DAT["Run"] == wildcards.run]
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
        fq = INDIR + "/{run}.fastq.gz",
        fa1 = rules.get_barcodes.output.fa1,
        fa2 = rules.get_barcodes.output.fa2
    output:
        tsv = OUTDIR + "/fbilr/{run}.tsv.gz"
    log:
        OUTDIR + "/fbilr/{run}.log"
    threads:
        THREADS
    shell:
        """(
        fbilr -t {threads} -w 200 -b {input.fa1},{input.fa2} {input.fq} \
            | pigz -p {threads} -c > {output} ) &> {log}
        """

rule report_matrix_summary:
    input:
        tsv = rules.fbilr.output.tsv
    output:
        tsv = OUTDIR + "/fbilr_stats/{run}.tsv"
    shell:
        """
        zcat {input} | awk '$2>=400' \
            | awk -v OFS=',' '{{print $3,$4,$5,$8,$9,$10,$11,$14}}' \
            | sort | uniq -c | awk -v OFS=',' '{{print $2,$1}}' \
            | sed 's/,/\\t/g' > {output}
        """

rule split_reads:
    input:
        fq = INDIR + "/{run}.fastq.gz",
        tsv1 = rules.fbilr.output,
        tsv2 = rules.get_barcodes.output.tsv
    output:
        directory(OUTDIR + "/splitted/{run}")
    log:
        OUTDIR + "/splitted/{run}.log"
    shell:
        """
        ./scripts/demux/split_reads.py -e 6 -l 400 {input} {output} &> {log}
        """

rule combine_reads:
    input:
        rules.split_reads.output
    output:
        fq = OUTDIR + "/combined/{run}/{cell}.fastq.gz"
    threads:
        THREADS
    shell:
        """
        fq1="{input}/fastqs/{wildcards.cell}_F.fastq"
        fq2="{input}/fastqs/{wildcards.cell}_R.fastq"
        ( cat $fq1; cat $fq2 | ./scripts/demux/reverse_fastq.py ) 、
            | pigz -p {threads} -c > {output}
        """

rule trim_reads:
    input:
        fq = OUTDIR + "/combined/{run}/{cell}.fastq.gz"
    output:
        directory(OUTDIR + "/trimmed/{run}/{cell}")
    log:
        OUTDIR + "/trimmed/{run}/{cell}.log"
    shell:
        """
        ./scripts/demux/trim_reads.py {input} {output} &> {log}
        """
