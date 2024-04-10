#!/usr/bin/env runsnakemake
import json
NAME = "HG001_Cell_350"
if "name" in config:
    NAME = config["name"]        
CONF = "data/%s.json" % NAME
d = json.load(open(CONF))
SPECIES = d["Species"]
CELLS = d["Cells"]
CHROMS = d["Chroms"]
STRAIN = d["Strain"]
ROOT_DIR = "results/%s" % NAME
THREADS = 12
print("Name:", NAME)
print("Species:", SPECIES)
print("Strain:", STRAIN)
print("Chroms:", len(CHROMS))
print("Cells:", len(CELLS))

def get_raw_bam(cell):
    return "../1_NanoStrandSeq/results/mapping/minimap2/%s/%s.bam" % (cell.split(".")[0], cell)

def get_bam(cell):
    return "../1_NanoStrandSeq/results/mapping/final/%s/%s.bam" % (cell.split(".")[0], cell)

if SPECIES == "Human":
    GENOME_FASTA = "/lustre/grp/tfclab/chenzg/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa"
    GENOME_SIZES = "/lustre/grp/tfclab/chenzg/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.sizes"
elif SPECIES == "Mouse":
    GENOME_FASTA = "/lustre/grp/tfclab/chenzg/species/mus_musculus/GRCm38.p6/GRCm38.canonical.genome.fa"
    GENOME_SIZES = "/lustre/grp/tfclab/chenzg/species/mus_musculus/GRCm38.p6/GRCm38.canonical.genome.sizes"
else:
    assert False

if STRAIN == "HG001":
    SNP_VCF = "../GRCh38-HG001-Variant-Calls/results/benchmark_autosomal_v4.2.1_chrx_v3.3.2.vcf.gz"
    SNP_BED = "../GRCh38-HG001-Variant-Calls/results/benchmark_autosomal_v4.2.1_chrx_v3.3.2.bed.gz"
elif STRAIN == "B6D2F1":
    SNP_VCF = "/lustre/grp/tfclab/chenzg/species/mus_musculus/GRCm38-C57-DBA-Variant-Calls/results/C57BL_6NJ_DBA_2J.mgp.v5.snps.dbSNP142.vcf.gz"
    SNP_BED = "../1_NanoStrandseq/data/Mouse_pseudo_high_confident_regions.bed.gz"
else:
    assert False
