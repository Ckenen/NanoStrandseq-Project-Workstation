#!/usr/bin/env python
import sys
import os
import gzip
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from pyBioInfo.IO.File import BedFile


def get_sequence(path, chrom):
    with pysam.FastaFile(path) as fasta:
        sequence = fasta.fetch(chrom)
    return sequence


def get_giab_high_confidence(path, chrom, length):
    hcrs = np.zeros(length, dtype=np.int)
    with BedFile(path) as f:
        for obj in f:
            if obj.chrom != chrom:
                continue
            for pos in range(obj.start, obj.end):
                hcrs[pos] = 1
    return hcrs


def get_parental_base(path, chrom, length):
    bases_pat = ["."] * length
    bases_mat = ["."] * length
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            row = line.strip("\n").split("\t")
            if row[6] != "PASS":
                continue
            if row[0] != chrom:
                continue
            pos = int(row[1]) - 1
            ref = row[3]
            alts = row[4].split(",")
            alleles = [ref] + alts
            attris = dict()
            for k, v in zip(row[8].split(":"), row[9].split(":")):
                attris[k] = v
            gt = attris["GT"]
            b1 = alleles[int(gt[0])]
            b2 = alleles[int(gt[2])]
            if gt[1] == "|":
                # pat
                if ref != b1:
                    m = max(len(ref), len(b1))
                    for n in range(m):
                        if n < len(ref):
                            if n < len(b1):
                                if b1[n] != ref[n]: 
                                    bases_pat[pos + n] = b1[n]
                            else:
                                bases_pat[pos + n] = "-"
                # mat
                if ref != b2:
                    m = max(len(ref), len(b2))
                    for n in range(m):
                        if n < len(ref):
                            if n < len(b2):
                                if b2[n] != ref[n]: 
                                    bases_mat[pos + n] = b2[n]
                            else:
                                bases_mat[pos + n] = "-"
            else:
                continue
                # assert b1 != b2
    return bases_pat, bases_mat
    


def load_base_matrix(path):
    with open(path) as f:
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
            base_counter = defaultdict(int)
            for item in row[8].split(","):
                k, v = item.split(":")
                v = int(v)
                base_counter[k] = v
            cell_base_counter_list = []
            for items in row[9].split(";"):
                d = defaultdict(int)
                for item in items.split(","):
                    k, v = item.split(":")
                    v = int(v)
                    d[k] = v
                cell_base_counter_list.append(d)
            yield pos, [base, cell_count, read_count, score, conf, base_counter, cell_base_counter_list]
        

def main():
    infile1 = None # genome.fasta
    infile2 = None # benchmark.bed
    infile3 = None # benchmark.phased.vcf
    infile4 = None # chr1.matrix
    chrom = None # chr1
    outdir = None
    infile1, infile2, infile3, infile4, chrom, outdir = sys.argv[1:]
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    sequence = get_sequence(infile1, chrom)
    hcrs = get_giab_high_confidence(infile2, chrom, len(sequence))
    bases_pat, bases_mat = get_parental_base(infile3, chrom, len(sequence))
    base_matrix_loader = load_base_matrix(infile4)
                
    # compare
    
    data1 = []
    for i in range(1, 6):
        counter1 = defaultdict(int) # pat
        counter2 = defaultdict(int) # mat
        parentals = []
        data1.append([counter1, counter2, parentals])
    cell_numbers = []
    
    current_pos_item = None
    eof = False
    
    for pos, (ref, hc, base_pat, base_mat) in enumerate(zip(sequence, hcrs, bases_pat, bases_mat)):
        # if pos < 10000000:
        #     continue
        # if pos % 1e6 == 0:
        #     print(pos)
    #     if pos >= 1e8:
    #         break
        if base_pat == ".":
            base_pat = ref
        if base_mat == ".":
            base_mat = ref
        base, cell_count, read_count, score, conf, base_counter, cell_base_counter_list = ".", 0, 0, 0, False, None, None
        
        item = None
        if current_pos_item is None:
            if eof:
                item = None
            else:
                try:
                    current_pos_item = next(base_matrix_loader)
                except StopIteration:
                    item = None
                    eof = True
        if current_pos_item is not None:
            if current_pos_item[0] > pos:
                item = None
            elif current_pos_item[0] == pos:
                item = current_pos_item[1]
                current_pos_item = None
            else:
                assert False
        
        if item:
            base, cell_count, read_count, score, conf, base_counter, cell_base_counter_list = item
            
        cell_count1 = 0
        if cell_count >= 1:
            score = 0
            if len(cell_base_counter_list) == 1:
                d = cell_base_counter_list[0]
                v1 = d[base]
                v2 = sum(d.values())
                if v2 >= 2:
                    cell_count1 += 1
                if v1 / v2 >= 0.8 and v2 >= 4:
                    score = 1
            else:
                for d in cell_base_counter_list:
                    v1 = d[base]
                    v2 = sum(d.values())
                    if v2 >= 2:
                        cell_count1 += 1
                    if v1 >= 2 and v1 / v2 >= 0.8:
                        score += 1
        
        for cutoff in [1, 2, 3, 4, 5]:
            # conf = score >= cutoff and score / cell_count >= 0.75
            conf = score >= cutoff and score / cell_count1 >= 0.75
            counter1, counter2, parentals = data1[cutoff - 1]
            s1 = "%s%s" % (ref, base_pat)
            s2 = "%s%s" % (ref, base_mat)
            if hc != 1:
                s1 = s1.lower()
                s2 = s2.lower()
            s3 = "%s%s" % (ref, base)
            if not conf:
                s3 = s3.lower()
            counter1[(s1, s3)] += 1
            counter2[(s2, s3)] += 1
            if hc == 1 and conf and base_pat != base_mat and base_pat != "-" and base_mat != "-":
                if base == base_pat:
                    parentals.append("P")
                elif base == base_mat:
                    parentals.append("M")
                else:
                    parentals.append("U")
            if cutoff == 2 and hc == 1 and not conf:
                cell_numbers.append(cell_count)
    
    for cutoff in [1, 2, 3, 4, 5]:   
        print("=" * 80)
        print("Score:", cutoff)
        for parental in [0, 1]:
            print("-" * 80)
            
            counter = data1[cutoff - 1][parental]
            parentals = data1[cutoff - 1][2]
            ks = []
            for k in counter.keys():
                ks.append(k[0])
                ks.append(k[1])
            ks = list(sorted(set(ks)))
            m = np.zeros((len(ks), len(ks)), dtype=np.int)
            for i in range(len(ks)):
                for j in range(len(ks)):
                    m[i][j] = counter[(ks[i], ks[j])]

            d = pd.DataFrame(m)
            d.columns = ks
            d.index = ks

            ks1 = list(filter(lambda item: item[0].isupper(), ks))
            ks2 = list(filter(lambda item: item[0].islower() and item[1] != ".", ks))
            ks3 = list(filter(lambda item: item[0].islower() and item[1] == ".", ks))
            ks4 = ks1 + ks2 + ks3
            d = d.loc[ks4][ks4]

            # pie
            if parental == 0:
                tmp1 = d.loc[ks1][ks1]
                tmp2 = d.loc[ks1][ks2]
                tmp3 = d.loc[ks1][ks3]
                tmp4 = d.loc[ks2][ks1]
                tmp5 = d.loc[ks2][ks2]
                tmp6 = d.loc[ks2][ks3]
                ds = [tmp1, tmp2, tmp3, tmp4, tmp5, tmp6]
                ss = [tmp.sum().sum() for tmp in ds]
                plt.figure(figsize=(4, 4))
                plt.pie(ss, labels=["H-H","H-L", "H-X", "L-H","L-L", "L-X"], autopct="%.1f%%")
                plt.tight_layout()
                plt.savefig(outdir + "/pie.cutoff%d.pdf" % cutoff, dpi=300)
                plt.close()
                
            # print("Precison of phasing:")
            if parental == 0:
                print("Precison of phasing:")
                counter = Counter(parentals)
                total = sum(counter.values())
                for k in ["P", "M", "U"]:
                    v = counter[k]
                    p = np.divide(v, total)
                    print(k, v, p, sep="\t")
            
            print("-" * 80)
            print("Parental: %d" % parental)
            # H-H
            print("Precision of assembly:")
            ks5 = list(filter(lambda item: "-" not in item, ks1))
            tmp = d.loc[ks5][ks5]
            v1, v2, v3 = 0, 0, 0
            for i in range(len(tmp)):
                for j in range(len(tmp.columns)):
                    s1 = tmp.index.values[i]
                    s2 = tmp.columns[j]
                    v = tmp.values[i][j]
                    if s1 == s2:
                        if s1[0] == s1[1]:
                            v1 += v
                        else:
                            v2 += v
                    else:
                        v3 += v
            print("Count1:", v1)
            print("Count2:", v2)
            print("Count3:", v3)
            p1 = (v1 + v2) / (v1 + v2 + v3)
            p2 = v2 / (v2 + v3)
            print("Precision:", p1)
            print("Precision:", p2)
            
        print("\n" * 5)
    

if __name__ == '__main__':
    main()
    