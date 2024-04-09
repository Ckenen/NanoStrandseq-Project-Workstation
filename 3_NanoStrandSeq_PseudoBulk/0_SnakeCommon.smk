#!/usr/bin/env runsnakemake

# NanoStrandSeq: 102.03x (63.65% duplicates)
# PacBio: 28.15x, 
# Ultralong: 34.30x

# NSS: downsample by NanoStrand-seq cells

NAMES = []
NAMES = ["PacBio.full", "Ultralong.full", "NSS.full"]
COVS = [2, 4, 5, 6, 8, 10, 12, 14, 15, 16, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]
# COVS = [90, 100]
REPEATS = [1, 2]

for source in ["PacBio", "Ultralong", "NSS", "NSS_RmDupFlag"]:
    for cov in COVS:
        if source == "Ultralong" and cov > 30:
            break
        if source == "PacBio" and cov > 25:
            break
        for rep in REPEATS:
            NAMES.append("%s.cov%d-r%d" % (source, cov, rep)) 
print("NAMES: %d" % len(NAMES))

THREADS = 12
OUTDIR = "results"
NSS_BAM_DIR = "data/nss_bams" # Directory of NanoStrand-seq bam files
FASTA = "/lustre/grp/tfclab/chenzg/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa"
SIZES = "/lustre/grp/tfclab/chenzg/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.sizes"
SNP_VCF = "/lustre/grp/tfclab/chenzg/repositories/GRCh38_HG001_SNP_Indel/GRCh38_HG001_SNP_Indel.GIAB_v4.2.1_and_v3.3.2.vcf.gz"
SNP_BED = "/lustre/grp/tfclab/chenzg/repositories/GRCh38_HG001_SNP_Indel/GRCh38_HG001_SNP_Indel.GIAB_v4.2.1_and_v3.3.2.bed.gz"
PBCCS_BAM = "../public/GIAB/HG001_GRCh38_NISTv4.2.1/PacBio_SequelII_CCS_11kb/HG001_GRCh38.haplotag.RTG.trio.bam"
ONTUL_BAM = "../public/GIAB/HG001_GRCh38_NISTv4.2.1/Ultralong_OxfordNanopore/NA12878-minion-ul_GRCh38.bam"
STRATIFICATION_TSV = "/lustre/grp/tfclab/chenzg/repositories/genome-stratifications/v3.1-genome-stratifications-GRCh38/v3.1-GRCh38-all-stratifications.tsv"
CLAIR_MODEL_CCS = "/lustre/grp/tfclab/chenzg/software/princess/bin/modules/ccs/model"
CLAIR_MODEL_ONT = "/lustre/grp/tfclab/chenzg/software/princess/bin/modules/ont/model"