#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('chaos', ["C1", "white", "C0"])

def main():
    group, outdir = sys.argv[1:]
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    infos = pd.read_excel("NanoStrandseq.xls")

    dat = infos[infos["Group"] == group]

    array = []
    for run, cell in dat[["Run", "Cell"]].values:
        path = "results/cwr/phased/%s/%s.txt" % (run, cell)
        d = pd.read_csv(path, sep="\t", index_col=0)
        s1 = d["Crick"].sum() 
        s2 = d["Watson"].sum()
        if s1 + s2 < 100000:
            continue
        s = d["Log2CWR"]
        s.name = cell
        array.append(s)

    dat = pd.DataFrame(array)
    dat = dat[dat.columns[:-3]]
    dat = dat[dat.abs().max(axis=1) >= 4]    
    matrix = dat
    values = matrix.abs().values.flatten()
    vmax = values[np.isfinite(values)].max()

    plt.figure(figsize=(10, 2 + len(matrix) * 0.2))
    plt.title("%s MAX(ABS(CWR))=%.2f" % (group, vmax))
    seaborn.heatmap(matrix, cmap=cmap, vmin=-vmax, vmax=vmax, annot=False, fmt=".1f")
    plt.tight_layout()
    plt.savefig(outdir + "/heatmap.png", dpi=300)
    plt.close()

    rows = []
    for chrom in dat.columns:
        dat1 = dat[dat[chrom].abs() < 0.5]
        rows.append([chrom, len(dat1), ",".join(dat1.index)])
        matrix = dat1
        values = matrix.abs().values.flatten()
        vmax = values[np.isfinite(values)].max()
        plt.figure(figsize=(10, 2 + len(matrix) * 0.2))
        plt.title("%s MAX(ABS(CWR))=%.2f" % (chrom, vmax))
        seaborn.heatmap(matrix, cmap=cmap, vmin=-vmax, vmax=vmax, annot=False, fmt=".1f")
        plt.tight_layout()
        plt.savefig(outdir + "/heatmap.%s.png" % chrom, dpi=300)
        plt.close()
    
    d = pd.DataFrame(rows)
    d.columns = ["Chrom", "Number", "Cells"]
    d.to_csv(outdir + "/stats.tsv", sep="\t", index=False)
    
    
if __name__ == '__main__':
    main()
    