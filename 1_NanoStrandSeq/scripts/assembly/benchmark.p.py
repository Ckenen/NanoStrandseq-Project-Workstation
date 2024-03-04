#!/usr/bin/env python
import sys
import os
import gzip
# import json
import multiprocessing
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt


def get_sequence(path, chrom, start, end):
    with pysam.FastaFile(path) as fasta:
        sequence = fasta.fetch(chrom, start, end)
    return sequence


def get_giab_high_confidence(path, chrom, start, end):
    length = end - start
    hcrs = np.zeros(length, dtype=np.int)
    with pysam.TabixFile(path) as f:
        try:
            for line in f.fetch(chrom, start, end):
                row = line.strip("\n").split("\t")
                s = int(row[1])
                e = int(row[2])
                for pos in range(max(s, start) - start, min(e, end) - start):
                    hcrs[pos] = 1
        except ValueError:
            pass
    return hcrs


def get_parental_base(path, chrom, start, end):
    length = end - start
    bases_pat = ["."] * length
    bases_mat = ["."] * length
    with pysam.VariantFile(path) as f:
        sample = list(f.header.samples)[0]
        for record in f.fetch(chrom, max(0, start - 50), end + 50):
            ps = record.samples[sample]["PS"]
            if ps != "PATMAT":
                continue
            ref = record.ref
            alleles = record.alleles
            gt = record.samples[sample]["GT"]
            pos = record.pos - 1
            idx = pos - start
            b1 = alleles[gt[0]]
            b2 = alleles[gt[1]]
            if b1 != ref: # pat
                m = max(len(ref), len(b1))
                for n in range(m):
                    if n < len(ref):
                        if n < len(b1):
                            if b1[n] != ref[n]:
                                idx0 = idx + n
                                if 0 <= idx0 < length: # indel
                                    bases_pat[idx0] = b1[n]
                        else:
                            idx0 = idx + n
                            if 0 <= idx0 < length:
                                bases_pat[idx0] = "-"
            if b2 != ref: # mat
                m = max(len(ref), len(b2))
                for n in range(m):
                    if n < len(ref):
                        if n < len(b2):
                            if b2[n] != ref[n]:
                                idx0 = idx + n
                                if 0 <= idx0 < length:
                                    bases_mat[idx0] = b2[n]
                        else:
                            idx0 = idx + n
                            if 0 <= idx0 < length:
                                bases_mat[idx0] = "-"
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
            while pos > require:
                yield null_item
                require += 1
            assert pos == require
            yield base, cell_count, read_count, score, conf, base_counter, cell_base_counter_list
            require += 1
    while require < end:
        yield null_item
        require += 1
    

def process(fasta, bed, vcf, mtx, chrom, start, end):
    # print(chrom, start, end, mtx, sep="\t")
    sequence = get_sequence(fasta, chrom, start, end)
    hcrs = get_giab_high_confidence(bed, chrom, start, end)
    bases_pat, bases_mat = get_parental_base(vcf, chrom, start, end)
    base_matrix_loader = load_base_matrix(mtx, chrom, start, end)
    
    data1 = []
    for i in range(1, 6): # min cell number
        counter1 = defaultdict(int) # pat
        counter2 = defaultdict(int) # mat
        parentals = []
        data1.append([counter1, counter2, parentals])
        
    cell_numbers = []
    
    for idx, (ref, hc, base_pat, base_mat, mtx_item) in enumerate(zip(sequence, hcrs, bases_pat, bases_mat, base_matrix_loader)):
        if base_pat == ".":
            base_pat = ref
        if base_mat == ".":
            base_mat = ref
        base, cell_count, read_count, score, conf, base_counter, cell_base_counter_list = mtx_item
                
        cell_count1 = 0
        if cell_count >= 1:
            score = 0
            if len(cell_base_counter_list) == 1:
                d = cell_base_counter_list[0]
                v1 = d[base]
                # print(base, v1)
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
    return data1


def main():
    fasta = None # genome.fasta
    bed = None # benchmark.bed
    vcf = None # benchmark.phased.vcf
    mtxdir = None # chr1.matrix
    chrom = None # chr1
    threads = 1 # threads
    outdir = None
    fasta, bed, vcf, mtxdir, chrom, threads, outdir = sys.argv[1:]
    threads = int(threads)
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    if threads == 1:
        array1 = []
        with open(os.path.join(mtxdir, "meta.tsv")) as f:
            for line in f:
                chrom, start, end, path = line.strip("\n").split("\t")
                fn = path.split("/")[-1]
                start = int(start)
                end = int(end)
                mtx = os.path.join(mtxdir, fn)
                array1.append(process(fasta, bed, vcf, mtx, chrom, start, end))
    else:
        pool = multiprocessing.Pool(threads)
        array = []
        with open(os.path.join(mtxdir, "meta.tsv")) as f:
            for line in f:
                chrom, start, end, path = line.strip("\n").split("\t")
                fn = path.split("/")[-1]
                start = int(start)
                end = int(end)
                mtx = os.path.join(mtxdir, fn)
                res = pool.apply_async(process, (fasta, bed, vcf, mtx, chrom, start, end))
                array.append(res)
        pool.close()
        pool.join()
        array1 = [item.get() for item in array]
        
    data1 = []
    for i in range(1, 6): # min cell number
        counter1 = defaultdict(int) # pat
        counter2 = defaultdict(int) # mat
        parentals = []
        
        for item in array1:
            for k, v in item[i - 1][0].items():
                counter1[k] += v
            for k, v in item[i - 1][1].items():
                counter2[k] += v
            for v in item[i - 1][2]:
                parentals.append(v)
        
        data1.append([counter1, counter2, parentals])
        
    # with open(outdir + "/details.json", "w+") as fw:
    #     json.dump(data1, fw, indent=4)
        
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
                # plt.show()
                plt.savefig(outdir + "/pie.cutoff%d.pdf" % cutoff, dpi=300)
                plt.close()

            # print("Precison of phasing:")
            if parental == 0:
                print("Precison of phasing:")
                counter = Counter(parentals)
                total = sum(counter.values())
                with open(outdir + "/cell.%d.txt" % cutoff, "w+") as fw:
                    for k in ["P", "M", "U"]:
                        v = counter[k]
                        p = np.divide(v, total)
                        print(k, v, p, sep="\t")
                        fw.write("\t".join(map(str, [k, v, p])) + "\n")
                    

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
    