#!/usr/bin/env runsnakemake

assembly_name = "HG001_Cell_350"
if "name" in config:
    assembly_name = config["name"]        
assembly_conf = "data/%s.json" % assembly_name
d = json.load(open(assembly_conf))
species = d["Species"]
cells = d["Cells"]
chroms = d["Chroms"]
assembly_dir = "results/%s" % assembly_name
if species == "Human":
    cellline = "HG001"
elif species == "Mouse":
    cellline = "Mouse"
else:
    assert False

print("Name:", assembly_name)
print("Species:", species)
print("Cell line:", cellline)
print("Chroms:", len(chroms))
print("Cells:", len(cells))

threads = 12

def get_raw_bam(cell):
    return "../1_NanoStrandSeq/results/mapping/minimap2/%s/%s.bam" % (cell.split(".")[0], cell)

def get_bam(cell):
    return "../1_NanoStrandSeq/results/mapping/mark_parental/%s/%s.bam" % (cell.split(".")[0], cell)

## Files

GENOMES = {
    "Human": {
        "GENOME_FASTA": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa",
        "GENOME_SIZE": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.sizes",
        "GENOME_MMI": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.mm2.map-ont.mmi",
        # "BLANK_BED": "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.Ns.bed.gz",
    },
    "Mouse": {
        "GENOME_FASTA": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.genome.fa",
        "GENOME_SIZE": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.genome.sizes",
        "GENOME_MMI": "/home/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.mm2.map-ont.mmi",
        # "BLANK_BED": "/date/chenzonggui/species/mus_musculus/GRCm38.p6/GRCm38.canonical.genome.Ns.bed.gz",
    }
}

BENCHMARKS = {
    "HG001": {
        "VCF": "../GRCh38-HG001-Variant-Calls/results/benchmark_autosomal_v4.2.1_chrx_v3.3.2.vcf.gz",
        "BED": "../GRCh38-HG001-Variant-Calls/results/benchmark_autosomal_v4.2.1_chrx_v3.3.2.bed.gz"
    },
    # "HG002": {
    #     "VCF": "../GIAB/HG002/HG002_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.corrected.vcf.gz",
    #     "BED": "../GIAB/HG002/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.sorted.bed.gz",
    # },
    "Mouse": {
        "VCF": "/home/chenzonggui/species/mus_musculus/GRCm38-C57-DBA-Variant-Calls/results/C57BL_6NJ_DBA_2J.mgp.v5.snps.dbSNP142.vcf.gz",
        "BED": "../1_NanoStrandseq/data/Mouse_pseudo_high_confident_regions.bed.gz"
    }
}
