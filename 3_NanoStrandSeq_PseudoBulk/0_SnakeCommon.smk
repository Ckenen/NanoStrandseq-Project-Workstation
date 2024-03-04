#!/usr/bin/env runsnakemake

# NanoStrandSeq: 102.03x (63.65% duplicates)
# PacBio: 28.15x, 
# Ultralong: 34.30x

# NSS # downsample by NanoStrand-seq cells

names = []
names = ["PacBio.full", "Ultralong.full", "NSS.full"]
covs = [2, 4, 5, 6, 8, 10, 12, 14, 15, 16, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]
# covs = [90, 100]
repeats = [1, 2]

for source in ["PacBio", "Ultralong", "NSS", "NSS_RmDupFlag"]:
    for cov in covs:
        if source == "Ultralong" and cov > 30:
            break
        if source == "PacBio" and cov > 25:
            break
        for rep in repeats:
            names.append("%s.cov%d-r%d" % (source, cov, rep))
            
print("Names: %d" % len(names))

threads = 12
THREADS = 12
outdir = "results"

# Directory of NanoStrand-seq bam files
NSS_BAM_DIR = "data/nss_bams"
FASTA = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa"
SIZES = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.sizes"
BENCHMARK_VCF = "../GRCh38-HG001-Variant-Calls/results/HG001_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.revised_mhc.vcf.gz"
BENCHMARK_BED = "../GRCh38-HG001-Variant-Calls/results/benchmark_autosomal_v4.2.1_chrx_v3.3.2.bed.gz"
