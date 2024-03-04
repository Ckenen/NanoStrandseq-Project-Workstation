#!/usr/bin/env
import pandas as pd

runs = [
    "20160928_PRJEB14185_CEU",
    "20210729_PRJNA742746_HG001"
]
dat = pd.read_excel("data/StrandSeq.xlsx")
dat["RunCell"] = ["%s/%s" % (run, cell) for run, cell in dat[["Run", "Cell"]].values]

# print(dat)

tmp = dat[[run in runs for run in dat["Run"]]]
run_cells = tmp["RunCell"]

def get_species(cell):
    return dat[dat["Cell"] == cell]["Species"].values[0]

def get_cellline(cell):
    return dat[dat["Cell"] == cell]["CellLine"].values[0]

GENOMES = {
    "Human": {
        "GENOME_FASTA": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa",
        "GENOME_SIZE": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.sizes",
        "GENOME_BWA": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.bwa.index",
    },
    "Mouse": {
        "GENOME_FASTA": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.genome.fa",
        "GENOME_SIZE": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.genome.sizes",
        "GENOME_BWA": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.bwa.index",
    }
}

BENCHMARKS = {
    "HG001": {
        "VCF": "../public/GRCh38-HG001-Variant-Calls/results/benchmark_autosomal_v4.2.1_chrx_v3.3.2.vcf.gz",
        "BED": "../public/GRCh38-HG001-Variant-Calls/results/benchmark_autosomal_v4.2.1_chrx_v3.3.2.bed.gz"
    },
    "HG002": {
        "VCF": "../public/GIAB/HG002/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.corrected.vcf.gz",
        "BED": "../public/GIAB/HG002/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.sorted.bed.gz",
    },
    "Mouse": {
        "VCF": "/home/chenzonggui/species/mus_musculus/mgp/C57BL_6NJ_DBA_2J.mgp.v5.snps.dbSNP142.vcf.gz",
        "BED": "../1_NanoStrandseq/data/Mouse_pseudo_high_confident_regions.bed.gz"
    }
}