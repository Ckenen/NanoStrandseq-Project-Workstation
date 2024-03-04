#!/usr/bin/env python
import sys
import sys
import os
import gzip
import multiprocessing
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import pysam
import matplotlib
# matplotlib.use("agg")
import matplotlib.pyplot as plt
from pyBioInfo.IO.File import BedFile
from pyBioInfo.Utils import ShiftLoader

def main():
    # infile1 = "../GIAB/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.sorted.bed.gz"
    # infile2 = "../GIAB/HG001/HG001_GRCh38.phased.patmat.bed.gz"
    # infile3 = "results/assembly/HG001_Cell50/inversions/inversions.bed.gz"
    # infile4 = "results/assembly/HG001_Cell50/round1/hets.all.bed.gz"
    infile1, infile2, infile3, infile4, outdir = sys.argv[1:]
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    regions = defaultdict(list)
    with BedFile(infile1) as f:
        for record in f:
            regions[record.chrom].append(record)
            
    pms = defaultdict(list)
    with BedFile(infile2) as f:
        for record in f:
            pms[record.chrom].append(record)
            
    for chrom, records in pms.items():
        loader = ShiftLoader(regions[chrom])
        records1 = []
        for record in records:
            if len(list(loader.fetch(obj=record))) > 0:
                records1.append(record)
        pms[chrom] = records1
        
    inversions = defaultdict(list)
    with BedFile(infile3) as f:
        for record in f:
            record.ratio = float(record.name.split(";")[2])
            inversions[record.chrom].append(record)
            
    hets = defaultdict(list)
    with BedFile(infile4) as f:
        for record in f:
            hets[record.chrom].append(record)
            
    # mark inversion
    for chrom, records in hets.items():
        loader = ShiftLoader(inversions[chrom])
        for record in records:
            ratio = None
            for inv in loader.fetch(obj=record):
                ratio = inv.ratio
            record.ratio = ratio
            
    # mark high confidence
    for chrom, records in hets.items():
        records1 = []
        loader = ShiftLoader(regions[chrom])
        for record in records:
            if len(list(loader.fetch(obj=record))) > 0:
                records1.append(record)
        hets[chrom] = records1
        
    results = dict()
    for chrom, records in hets.items():
        n = len(pms[chrom])
        matrix = np.zeros((4, 3), dtype=np.int)
        loader = ShiftLoader(pms[chrom])
        for obj1 in records:
            obj2 = list(loader.fetch(obj=obj1))
            g = 1
            if obj1.ratio is None:
                g = 1
            else:
                if obj1.ratio >= 0.9:
                    g = 2
                else:
                    g = 3
            t = 2
            if len(obj2) == 1:
                obj2 = obj2[0]
                if obj1.name[0] == obj2.name[0] and obj1.name[2] == obj2.name[2]:
                    t = 0
                elif obj1.name[0] == obj2.name[2] and obj1.name[2] == obj2.name[0]:
                    t = 1
            matrix[0][t] += 1
            matrix[g][t] += 1
        d = pd.DataFrame(matrix)
        d.columns = ["HP1.Pat", "HP1.Mat", "Other"]
        d.index = ["All", "Not.Inv", "Hom.Inv", "Het.Inv"]
        d["GIAB"] = n
        d["Precision"] = d[["HP1.Pat", "HP1.Mat"]].max(axis=1) / d[["HP1.Pat", "HP1.Mat"]].sum(axis=1)
        d.to_csv(os.path.join(outdir, "%s.tsv" % chrom), sep="\t")
        results[chrom] = d
        
    rows = []
    for chrom, d in results.items():
        vs = d.values
        v1 = vs[1][0] + vs[2][1]
        v2 = vs[1][1] + vs[2][0]
        v3 = vs[0][3]
        p = max(v1, v2) / (v1 + v2)
        recall = max(v1, v2) / v3
        rows.append([chrom, v1, v2, v3, p, recall])
    d = pd.DataFrame(rows)
    d.columns = ["Chrom", "HP1.Pat", "HP1.Mat", "GIAB", "Precision", "Recall"]
    d.to_csv(os.path.join(outdir, "overall.tsv"), sep="\t", index=False)
            
if __name__ == '__main__':
    main()
    