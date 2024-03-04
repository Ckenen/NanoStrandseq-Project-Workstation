#!/usr/bin/env python
import os
import sys
import gzip
from collections import defaultdict
import numpy as np
import pysam
import multiprocessing
from pyBioInfo.IO.File import BedFile

def get_sequence(path, chrom, start, end):
    with pysam.FastaFile(path) as fasta:
        sequence = fasta.fetch(chrom, start, end)
    return sequence


def get_giab_high_confidence(path, chrom, start, end):
    length = end - start
    hcrs = np.zeros(length, dtype=np.int)
    # with pysam.TabixFile(path) as f:
    #     try:
    #         for line in f.fetch(chrom, start, end):
    #             row = line.strip("\n").split("\t")
    #             s = int(row[1])
    #             e = int(row[2])
    #             for pos in range(max(s, start) - start, min(e, end) - start):
    #                 hcrs[pos] = 1
    #     except ValueError:
    #         pass
    return hcrs


def get_parental_base(path, chrom, start, end):
    length = end - start
    bases_pat = ["."] * length
    bases_mat = ["."] * length
    # with pysam.VariantFile(path) as f:
    #     sample = list(f.header.samples)[0]
    #     for record in f.fetch(chrom, max(0, start - 50), end + 50):
    #         ps = record.samples[sample]["PS"]
    #         if ps != "PATMAT":
    #             continue
    #         ref = record.ref
    #         alleles = record.alleles
    #         gt = record.samples[sample]["GT"]
    #         pos = record.pos - 1
    #         idx = pos - start
    #         b1 = alleles[gt[0]]
    #         b2 = alleles[gt[1]]
    #         if b1 != ref: # pat
    #             m = max(len(ref), len(b1))
    #             for n in range(m):
    #                 if n < len(ref):
    #                     if n < len(b1):
    #                         if b1[n] != ref[n]:
    #                             idx0 = idx + n
    #                             if 0 <= idx0 < length: # indel
    #                                 bases_pat[idx0] = b1[n]
    #                     else:
    #                         idx0 = idx + n
    #                         if 0 <= idx0 < length:
    #                             bases_pat[idx0] = "-"
    #         if b2 != ref: # mat
    #             m = max(len(ref), len(b2))
    #             for n in range(m):
    #                 if n < len(ref):
    #                     if n < len(b2):
    #                         if b2[n] != ref[n]:
    #                             idx0 = idx + n
    #                             if 0 <= idx0 < length:
    #                                 bases_mat[idx0] = b2[n]
    #                     else:
    #                         idx0 = idx + n
    #                         if 0 <= idx0 < length:
    #                             bases_mat[idx0] = "-"
    return bases_pat, bases_mat


def load_base_matrix(path, chrom, start, end):
    require = start
    null_item = ".", 0, 0, 0, False, None, None
    with gzip.open(path, "rt") as f:
        for line in f:
            row = line.strip("\n").split("\t")
            pos = int(row[0])
            ref_base = row[1]
            base = row[2]
            base_read_count = int(row[3])
            cell_count = int(row[4])
            read_count = int(row[5])
            score = int(row[6])
            conf = row[7] == "T"
#             base_counter = defaultdict(int)
#             for item in row[8].split(","):
#                 k, v = item.split(":")
#                 v = int(v)
#                 base_counter[k] = v
#             cell_base_counter_list = []
#             for items in row[9].split(";"):
#                 d = defaultdict(int)
#                 for item in items.split(","):
#                     k, v = item.split(":")
#                     v = int(v)
#                     d[k] = v
#                 cell_base_counter_list.append(d)
            while pos > require:
                yield null_item
                require += 1
            assert pos == require
            yield base, cell_count, read_count, score, conf, row[8], row[9]
            require += 1
    while require < end:
        yield null_item
        require += 1

def process(fasta, bed, vcf, mtx1, mtx2, chrom, start, end, outfile):
    # print(chrom, start, end, outfile)
    sequence = get_sequence(fasta, chrom, start, end)
    hcrs = get_giab_high_confidence(bed, chrom, start, end)
    bases_pat, bases_mat = get_parental_base(vcf, chrom, start, end)    
    loader1 = load_base_matrix(mtx1, chrom, start, end)
    loader2 = load_base_matrix(mtx2, chrom, start, end)
    
    fw = gzip.open(outfile, "wt")
    for idx, (ref, hc, base_pat, base_mat, item1, item2) in enumerate(zip(sequence, hcrs, bases_pat, bases_mat, loader1, loader2)):
        pos = idx + start
        base1, cell_count1, read_count1, score1, conf1, base_counter1, cell_base_counter_list1 = item1
        base2, cell_count2, read_count2, score2, conf2, base_counter2, cell_base_counter_list2 = item2
        s1 = "."
        s2 = "."
        if cell_base_counter_list1 is not None:
            s1 = cell_base_counter_list1
        if cell_base_counter_list2 is not None:
            s2 = cell_base_counter_list2
        # print(pos, ref, base_pat, base_mat, hc, base1, base2, s1, s2, sep="\t")
        line = "\t".join(map(str, [pos, ref, base_pat, base_mat, hc, base1, base2, s1, s2]))
        fw.write(line + "\n")
    fw.close()
    

def main():    
    # fasta = "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa"
    # bed = "../GIAB/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.sorted.bed.gz"
    # vcf = "../GIAB/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark_hifiasm_v11_phasetransfer.corrected.vcf.gz"
    # mtxdir1 = "results/assembly/HG001_Cell50/round1/matrix/chr22.hp1"
    # mtxdir2 = "results/assembly/HG001_Cell50/round1/matrix/chr22.hp2"
    # chrom = "chr22"
    # threads = 10
    # outdir = "merged_out"
    fasta, bed, vcf, mtxdir1, mtxdir2, chrom, threads, outdir = sys.argv[1:]
    threads = int(threads)
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    pool = None
    if threads > 1:
        pool = multiprocessing.Pool(threads)
    results = []
    
    with open(os.path.join(mtxdir1, "meta.tsv")) as f, open(os.path.join(outdir, "meta.tsv"), "w+") as fw:
        for line in f:
            chrom, start, end, path = line.strip("\n").split("\t")
            start = int(start)
            end = int(end)
            fname = path.split("/")[-1]
            fw.write("\t".join(map(str, [chrom, start, end, fname])) + "\n")
            mtx1 = os.path.join(mtxdir1, fname)
            mtx2 = os.path.join(mtxdir2, fname)
            outfile = os.path.join(outdir, fname)
            if pool is None:
                process(fasta, bed, vcf, mtx1, mtx2, chrom, start, end, outfile)
            else:
                res = pool.apply_async(process, (fasta, bed, vcf, mtx1, mtx2, chrom, start, end, outfile))
                results.append(res)
    
    if pool is not None:
        pool.close()
        pool.join()
        for res in results:
            assert res.successful()


if __name__ == '__main__':
    main()
    