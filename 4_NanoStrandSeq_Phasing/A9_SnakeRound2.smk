#!/usr/bin/env runsnakemake
include: "0_SnakeNanoStrandSeq.smk"
outdir = assembly_dir + "/round2"
hps = ["hp1", "hp2"]
# chroms = ["chr22"]

rule all:
    input:
        #expand(outdir + "/splited_haplotype/{cell}", cell=cells),
        #expand(outdir + "/bams/{chrom}.{hp}.bam", chrom=chroms, hp=hps),
        #outdir + "/merged.bam",
        #expand(outdir + "/matrix/{chrom}.{hp}", chrom=chroms, hp=hps),
        expand(outdir + "/benchmark/{chrom}.{hp}", chrom=chroms, hp=hps),
        #expand(outdir + "/matrix2/{chrom}", chrom=chroms),
        #expand(outdir + "/matrix2_reduced/{chrom}", chrom=chroms),
        #expand(outdir + "/hets/{chrom}.bed", chrom=chroms),
        #outdir + "/hets.all.bed.gz",
        outdir + "/hets.benchmark",
        #expand(outdir + "/cuteSV/{chrom}.{hp}.vcf.gz", chrom=chroms, hp=hps),
        expand(outdir + "/cuteSV.{hp}.vcf.gz", hp=hps),
        expand(outdir + "/cutesv/quantify/{chrom}.{hp}.on_{hpA}.txt.gz", chrom=chroms, hp=hps, hpA=hps),
        expand(outdir + "/cutesv/quantify_lite/{chrom}.{hp}.on_{hpA}.tsv", chrom=chroms, hp=hps, hpA=hps),

def get_bam(cell):
    run = cell.split(".")[0]
    return "results/mapping/final/%s/%s.bam" % (run, cell)

rule split_haplotype:
    input:
        bam = lambda wildcards: get_bam(wildcards.cell),
        bed = assembly_dir + "/round1/hets.all.bed.gz"
    output:
        out = directory(outdir + "/splited_haplotype/{cell}")
    log:
        outdir + "/splited_haplotype/{cell}.log"
    threads:
        4
    shell:
        """
        ./scripts/assembly/split_haplotype.v1.py {input.bam} {input.bed} {output.out} &> {log}
        """

rule make_haplotype_bam:
    input:
        expand(rules.split_haplotype.output.out, cell=cells),
    output:
        tmp1 = temp(outdir + "/bams/{chrom}.hp1.tmp.bam"),
        tmp2 = temp(outdir + "/bams/{chrom}.hp2.tmp.bam"),
        bam1 = outdir + "/bams/{chrom}.hp1.bam",
        bam2 = outdir + "/bams/{chrom}.hp2.bam",
        bai1 = outdir + "/bams/{chrom}.hp1.bam.bai",
        bai2 = outdir + "/bams/{chrom}.hp2.bam.bai"
    log:
        outdir + "/bams/{chrom}.log"
    threads:
        4
    shell:
        """(
        ./scripts/strandtools/make_haplotype_bam_r2.py \
            {assembly_conf} {wildcards.chrom} {output.tmp1} {output.tmp2}
        samtools sort -@ {threads} -o {output.bam1} {output.tmp1}
        samtools sort -@ {threads} -o {output.bam2} {output.tmp2}
        samtools index -@ {threads} {output.bam1}
        samtools index -@ {threads} {output.bam2} ) &> {log}
        """

rule merge_bams:
    input:
        bams1 = [outdir + "/bams/%s.hp1.bam" % c for c in chroms],
        bams2 = [outdir + "/bams/%s.hp2.bam" % c for c in chroms],
    output:
        bam = outdir + "/merged.bam"
    threads:
        12
    shell:
        """
        (   
            samtools view -H --no-PG {input.bams1[0]} | grep -v '^@PG'
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

rule generate_base_matrix:
    input:
        bam = outdir + "/bams/{chrom}.{hp}.bam",
        fa = lambda wildcards: get_genome_fasta(cells[0]),
    output:
        mtx = directory(outdir + "/matrix/{chrom}.{hp}")
    log:
        outdir + "/matrix/{chrom}.{hp}.log"
    threads:
        24
    shell:
        """
        ./scripts/assembly/generate_base_matrix.p.v1.py \
            {input.bam} {input.fa} {wildcards.chrom} {threads} {output.mtx} &> {log}
        """

rule benchmark:
    input:
        fa = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa",
        bed = "../GIAB/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.sorted.bed.gz",
        vcf = "../GIAB/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.corrected.vcf.gz",
        mtx = rules.generate_base_matrix.output.mtx
    output:
        out = directory(outdir + "/benchmark/{chrom}.{hp}")
    log:
        outdir + "/benchmark/{chrom}.{hp}.log"
    threads:
        24
    shell:
        """
        ./scripts/assembly/benchmark.p.py \
            {input.fa} {input.bed} {input.vcf} {input.mtx} \
            {wildcards.chrom} {threads} {output.out} &> {log}
        """

rule merge_base_matrix:
    input:
        fa = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa",
        bed = "../GIAB/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.sorted.bed.gz",
        vcf = "../GIAB/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.corrected.vcf.gz",
        mtxdir1 = outdir + "/matrix/{chrom}.hp1",
        mtxdir2 = outdir + "/matrix/{chrom}.hp2"
    output:
        directory(outdir + "/matrix2/{chrom}")
    log:
        log = outdir + "/matrix2/{chrom}.log"
    threads:
        12
    shell:
        """
        ./scripts/assembly/merge_base_matrix.p.py {input.fa} {input.bed} {input.vcf} \
            {input.mtxdir1} {input.mtxdir2} {wildcards.chrom} {threads} {output} &> {log}
        """

rule reduce_matrix:
    input:
        outdir + "/matrix2/{chrom}"
    output:
        directory(outdir + "/matrix2_reduced/{chrom}")
    log:
        outdir + "/matrix2_reduced/{chrom}.log"
    threads:
        24
    shell:
        """
        ./scripts/strandtools/reduce_matrix2.py {input} {threads} {output} &> {log}
        """

# rule get_lite_matrix2:
#     input:
#         outdir + "/matrix2/{chrom}"
#     output:
#         directory(outdir + "/matrix2_lite/{chrom}")
#     threads:
#         4
#     shell:
#         """
#         mkdir {output}
#         for path in {input}/*.matrix.gz; do
#             bn=`basename $path`
#             pigz -p {threads} -d -c $path \
#                 | ./scripts/strandtools/get_lite_matrix2.py \
#                 | pigz -p {threads} -c > {output}/$bn
#         done
        # """

rule get_hets:
    input:
        bed = assembly_dir + "/snv/concat/nanocaller.vcf.gz",
        mtxdir = outdir + "/matrix2/{chrom}",
    output:
        bed = outdir + "/hets/{chrom}.bed",
    log:
        outdir + "/hets/{chrom}.log"
    threads:
        12
    shell:
        """
        ./scripts/assembly/get_hets.p.py {input.bed} {input.mtxdir} {wildcards.chrom} {threads} {output.bed} &> {log}
        """

rule merge_hets:
    input:
        beds = expand(rules.get_hets.output.bed, chrom=chroms)
    output:
        bed = outdir + "/hets.all.bed.gz",
        tbi = outdir + "/hets.all.bed.gz.tbi"
    shell:
        """
        cat {input.beds} | sort -k1,1 -k2,2n | bgzip -c > {output.bed}
        tabix -p bed {output.bed}
        """
        
rule benchmark_hets:
    input:
        infile1 = "../GIAB/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.sorted.bed.gz",
        infile2 = "../GIAB/HG001/HG001_GRCh38.phased.patmat.bed.gz",
        infile3 = assembly_dir + "/inversions/inversions.bed.gz",
        infile4 = rules.merge_hets.output.bed,
    output:
        directory(outdir + "/hets.benchmark")
    log:
        outdir + "/hets.benchmark.log"
    shell:
        """
        ./scripts/assembly/benchmark_hets.py {input} {output} &> {log}
        """

rule cutesv:
    input:
        bam = outdir + "/bams/{chrom}.{hp}.bam",
        bai = outdir + "/bams/{chrom}.{hp}.bam.bai",
        fa = lambda wildcards: get_genome_fasta(cells[0])
    output:
        wd = temp(directory(outdir + "/cuteSV/{chrom}.{hp}.wd")),
        vcf = temp(outdir + "/cuteSV/{chrom}.{hp}.vcf"),
        vcf2 = outdir + "/cuteSV/{chrom}.{hp}.vcf.gz",
        tbi2 = outdir + "/cuteSV/{chrom}.{hp}.vcf.gz.tbi"
    log:
        outdir + "/cuteSV/{chrom}.{hp}.log"
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
        vcfs = [outdir + "/cuteSV/%s.{hp}.vcf.gz" % c for c in chroms]
    output:
        vcf = outdir + "/cuteSV.{hp}.vcf.gz",
        tbi = outdir + "/cuteSV.{hp}.vcf.gz.tbi"
    shell:
        """
        bcftools concat -a {input.vcfs} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule quantify_sv:
    input:
        vcf = outdir + "/cuteSV/{chrom}.{hp}.vcf.gz",
        bam = outdir + "/bams/{chrom}.{hpA}.bam"
    output:
        tmp = temp(outdir + "/cutesv/quantify/{chrom}.{hp}.on_{hpA}.txt"),
        txt = outdir + "/cutesv/quantify/{chrom}.{hp}.on_{hpA}.txt.gz"
    threads:
        threads
    shell:
        """
        ../6_SNV_SV_Comparison/scripts/quantify_sv.p.py {input.vcf} {threads} {input.bam} {output.tmp}
        awk 'NR==1||$1=="{wildcards.chrom}"' {output.tmp} | pigz -p {threads} -c > {output.txt}
        """

rule lite_quantify:
    input:
        txt = outdir + "/cutesv/quantify/{source}.txt.gz"
    output:
        txt = outdir + "/cutesv/quantify_lite/{source}.tsv"
    shell:
        """
        gzip -d -c {input.txt} \
            | awk -v OFS='\\t' '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}' > {output.txt}
        """

# rule stat_bin_read_count:
#     input:
#         bam = outdir + "/mark_haplotype/{run}/{cell}.bam"
#         # bai = outdir + "/mark_haplotype/{run}/{cell}.bam.bai"
#     output:
#         out = directory(outdir + "/barplot/{run}/{cell}")
#     log:
#         outdir + "/barplot/{run}/{cell}.log"
#     shell:
#         """
#         ./scripts/strandtools/stat_bin_read_count.py {input.bam} {output.out} &> {log}
#         """

# def get_pdfs(wildcards):
#     run = wildcards.run
#     cells = dat[dat["Run"] == run]["Cell"]
#     paths = []
#     for cell in cells:
#         path = outdir + "/barplot/%s/%s" % (run, cell)
#         paths.append(path)
#     return paths

# rule merge_pdf:
#     input:
#         lambda wildcards: get_pdfs(wildcards)
#     output:
#         pdf1 = outdir + "/barplot/{run}.all_cells.pdf",
#         pdf2 = outdir + "/barplot/{run}.all_cells.Trim.pdf",
#         pdf3 = outdir + "/barplot/{run}.all_cells.IsHC.pdf",
#         pdf4 = outdir + "/barplot/{run}.all_cells.IsHC.Trim.pdf",
#         pdf5 = outdir + "/barplot/{run}.all_cells.RmDup.pdf",
#         pdf6 = outdir + "/barplot/{run}.all_cells.RmDup.Trim.pdf",
#         pdf7 = outdir + "/barplot/{run}.all_cells.RmDup.IsHC.pdf",
#         pdf8 = outdir + "/barplot/{run}.all_cells.RmDup.IsHC.Trim.pdf"
#     params:
#         prefix = outdir + "/barplot/{run}.all_cells"
#     shell:
#         """
#         params=""; for d in {input}; do params="$params $d/barplot.pdf"; done
#         merge_pdf.py $params {params.prefix}.pdf
#         params=""; for d in {input}; do params="$params $d/barplot.Trim.pdf"; done
#         merge_pdf.py $params {params.prefix}.Trim.pdf
#         params=""; for d in {input}; do params="$params $d/barplot.IsHC.pdf"; done
#         merge_pdf.py $params {params.prefix}.IsHC.pdf
#         params=""; for d in {input}; do params="$params $d/barplot.IsHC.Trim.pdf"; done
#         merge_pdf.py $params {params.prefix}.IsHC.Trim.pdf
#         params=""; for d in {input}; do params="$params $d/barplot.RmDup.pdf"; done
#         merge_pdf.py $params {params.prefix}.RmDup.pdf
#         params=""; for d in {input}; do params="$params $d/barplot.RmDup.Trim.pdf"; done
#         merge_pdf.py $params {params.prefix}.RmDup.Trim.pdf
#         params=""; for d in {input}; do params="$params $d/barplot.RmDup.IsHC.pdf"; done
#         merge_pdf.py $params {params.prefix}.RmDup.IsHC.pdf
#         params=""; for d in {input}; do params="$params $d/barplot.RmDup.IsHC.Trim.pdf"; done
#         merge_pdf.py $params {params.prefix}.RmDup.IsHC.Trim.pdf
#         """

# # common rules

# rule sort_bam:
#     input:
#         bam = "{prefix}.bam"
#     output:
#         bam = "{prefix}.sorted.bam"
#     shell:
#         """
#         samtools sort -o {output.bam} {input.bam}
#         """

# rule bam_index:
#     input:
#         bam = "{prefix}.bam"
#     output:
#         bai = "{prefix}.bam.bai"
#     shell:
#         """
#         samtools index {input.bam}
#         """

# rule bam_flagstat:
#     input:
#         bam = "{prefix}.bam"
#     output:
#         txt = "{prefix}.flagstat"
#     shell:
#         """
#         samtools flagstat {input.bam} > {output.txt}
#         """
