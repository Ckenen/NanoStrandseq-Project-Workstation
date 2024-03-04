#!/usr/bin/env runsnakemake

genome_fasta = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa"

hifi_fastq_dir = "../public/PRJNA540705_HG001_PacBio_CCS"
hifi_samples = ["SRR9001768", "SRR9001769", "SRR9001770", "SRR9001771", "SRR9001772", "SRR9001773"]
hifi_fastqs = [hifi_fastq_dir + "/%s.fastq.gz" % s for s in hifi_samples]

def get_ccs_fastq(sample):
    return hifi_fastq_dir + "/%s.fastq.gz" % sample

contig_fasta = "results/assembly/wtdbg2/asm.cns.fa"

cluster_fasta = "results/ss2contig/tgs/cluster/clusters.fasta"
cluster_sizes = "results/ss2contig/tgs/cluster/clusters.sizes"
clusters = []
if os.path.exists(cluster_fasta + ".fai"):
    for line in open(cluster_fasta + ".fai"):
        row = line.strip("\n").split("\t")
        if int(row[1]) > 20000000:
            clusters.append(row[0])
# print(len(clusters))
# print(clusters)

import glob
strand_seq_cells = []
for path in glob.glob("../2_Strandseq/results/prepare/cutadapt/20160928_PRJEB14185_CEU/*_1.fastq.gz"):
    cell = path.split("/")[-1][:-11]
    strand_seq_cells.append(cell)

import json
nano_strand_seq_cells = json.load(open("../4_NanoStrandseq_Assembly/data/HG001_Cell_350.json"))["Cells"]
nano_strand_seq_cells.sort()
nano_strand_seq_cells = nano_strand_seq_cells[:132]

def get_strand_seq_fastqs(cell):
    path1 = "../2_Strandseq/results/prepare/cutadapt/20160928_PRJEB14185_CEU/%s_1.fastq.gz" % cell
    path2 = "../2_Strandseq/results/prepare/cutadapt/20160928_PRJEB14185_CEU/%s_2.fastq.gz" % cell
    return [path1, path2]

def get_nano_strand_seq_fastq(cell):
    path = "../1_NanoStrandseq/results/demux/trimmed/%s/%s.fastq.gz" % (cell.split(".")[0], cell)
    return path