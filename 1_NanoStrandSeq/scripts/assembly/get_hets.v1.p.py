#!/usr/bin/env python
import sys
import gzip
import os
import multiprocessing
from collections import defaultdict
import numpy as np
import pandas as pd
import pysam


def decode(s):
    ds = []
    if s == ".":
        return ds
    for items in s.split(";"):
        d = defaultdict(int)
        for item in items.split(","):
            k, v = item.split(":")
            v = int(v)
            d[k] = v
        ds.append(d)
    return ds


def get_score(base, ds):
    cell_count = 0
    score = 0
    if len(ds) >= 1:
        score = 0
        if len(ds) == 1:
            d = ds[0]
            v1 = d[base]
            v2 = sum(d.values())
            if v2 >= 2:
                cell_count += 1
            if v1 / v2 >= 0.8 and v2 >= 4:
                score = 1
        else:
            for d in ds:
                v1 = d[base]
                v2 = sum(d.values())
                if v2 >= 2:
                    cell_count += 1
                if v1 >= 2 and v1 / v2 >= 0.8:
                    score += 1
    return cell_count, score


def process(path, chrom, start, end, hets):   
    matrix = np.zeros((4, 2), dtype=np.int)
    phased_hets = []
    
    with gzip.open(path, "rt") as f:
        for line in f:
            row = line.strip("\n").split("\t")
            if row[7] == "." and row[8] == ".":
                continue
            pos = int(row[0])
            base_ref = row[1]
            base_pat = row[2] if row[2] != "." else base_ref
            base_mat = row[3] if row[3] != "." else base_ref
            hc = row[4] == "1"
            base_hp1, base_hp2 = row[5], row[6]
            ds1, ds2 = decode(row[7]), decode(row[8])
            cell_count1, score1 = get_score(base_hp1, ds1)
            cell_count2, score2 = get_score(base_hp2, ds2)
            conf1 = score1 >= 1 and score1 / cell_count1 >= 0.75
            conf2 = score2 >= 1 and score2 / cell_count2 >= 0.75
            
            keep = False
            if conf1:
                if conf2:
                    if base_hp1 != base_hp2:
                        keep = True
                else:
                    ret = hets.get(pos)
                    if ret is not None:
                        if base_hp1 == ret[0]:
                            base_hp2 = ret[1]
                            keep = True
                        elif base_hp1 == ret[1]:
                            base_hp2 = ret[0]
                            keep = True
            else:
                if conf2:
                    ret = hets.get(pos)
                    if ret is not None:
                        if base_hp2 == ret[0]:
                            base_hp1 = ret[1]
                            keep = True
                        elif base_hp2 == ret[1]:
                            base_hp1 = ret[0]
                            keep = True
            
            if keep and base_hp1 != "-" and base_hp2 != "-":
                phased_hets.append([chrom, pos, pos + 1, "%s|%s" % (base_hp1, base_hp2)])
                if hc:
                    if base_pat == base_hp1:
                        matrix[0][0] += 1
                    else:
                        matrix[0][1] += 1
                    if base_pat == base_hp2:
                        matrix[1][0] += 1
                    else:
                        matrix[1][1] += 1
                    if base_mat == base_hp1:
                        matrix[2][0] += 1
                    else:
                        matrix[2][1] += 1
                    if base_mat == base_hp2:
                        matrix[3][0] += 1
                    else:
                        matrix[3][1] += 1
                        
    d = pd.DataFrame(matrix)
    d.columns = ["Same", "Diff"]
    d.index = ["Pat_HP1", "Pat_HP2", "Mat_HP1", "Mat_HP2"]
    d["Precision"] = d["Same"] / d.sum(axis=1)
    d.index.name = "Comparison"
    
    return phased_hets, d


def load_hets(path, chrom):
    hets = dict()
    with pysam.VariantFile(path) as f:
        name = list(f.header.samples)[0]
        for record in f.fetch(chrom):
            alleles = record.alleles
            gt = record.samples[name]["GT"]
            a1 = alleles[gt[0]]
            a2 = alleles[gt[1]]
            if a1 == a2 or len(a1) > 1 or len(a2) > 1:
                continue
            hets[record.pos - 1] = [a1, a2]
    return hets
    
    # hets = dict()
    # with pysam.TabixFile(path) as f:
    #     for line in f.fetch(chrom):
    #         row = line.strip("\n").split("\t")
    #         chrom, start, end, name = row[:4]
    #         base1, base2 = name[0], name[2]
    #         hets[int(start)] = [base1, base2]
    # return hets

def main():
    # hets.bed.gz mtxdir chrom 8 outfile
    bed, mtxdir, chrom, threads, outfile = sys.argv[1:]
    threads = int(threads)
    
    hets = load_hets(bed, chrom)
    
    results = []
    pool = multiprocessing.Pool(threads)
    with open(os.path.join(mtxdir, "meta.tsv")) as f:
        for line in f:
            chrom, start, end, fname = line.strip("\n").split("\t")
            start, end = int(start), int(end)
            mtx = os.path.join(mtxdir, fname)
            r = pool.apply_async(process, (mtx, chrom, start, end, hets))
            results.append(r)
    pool.close()
    pool.join()
    results = [r.get() for r in results]
        
    stats_all = None
    with open(outfile, "w+") as fw:
        for phased_hets, stats in results:
            if stats_all is None:
                stats_all = stats
            else:
                stats_all = stats_all + stats
            for row in phased_hets:
                fw.write("\t".join(map(str, row)) + "\n")
    stats_all["Precision"] = stats_all["Same"] / stats_all.sum(axis=1)
    stats_all.to_csv(sys.stdout, sep="\t")
    
    
if __name__ == '__main__':
    main()
    