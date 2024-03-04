#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt


def get_sharp(vs, cutoff1, cutoff2):
    ps = []
    for i in range(0, len(vs)):
        v = vs[i]
        if v < cutoff1:
            continue
        i1 = max(i - 5, 0)
        i2 = min(i + 6, len(vs))
        vs0 = vs[i1:i2]
        mean = np.mean(vs0)
        std = np.std(vs0)
        if std > 0:
            z = (v - mean) / std
            r = sum(vs0 < cutoff2) / len(vs0)
            if z > 2 and r > 0.5:
                ps.append(i)
    return ps


def main():
    infile, outdir = sys.argv[1:]
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    dat = pd.read_csv(infile, sep="\t")

    vs = dat["Crick"] + dat["Watson"]
    vs = vs[vs > 0]
    vs = np.sort(vs)

    c1 = vs[int(len(vs) * 0.1)]
    c2 = vs[int(len(vs) * 0.9)]
        
    if len(vs) > 10:
        xs = np.array(vs)
        ys = np.arange(len(xs))

        plt.figure(figsize=(4, 4))
        plt.plot(xs, ys)
        plt.axvline(c1, lw=1, ls="--", color="grey")
        plt.axvline(c2, lw=1, ls="--", color="grey")
        plt.tight_layout()
        plt.savefig(outdir + "/cumulative.png", dpi=300)
        plt.close()

    vs1 = vs[(vs > c1) & (vs < c2)]
    if len(vs1) == 0:
        vs1 = vs
    mean = np.mean(vs1)
    std = np.std(vs1)
    cutoff1 = mean * 0.3
    cutoff2 = mean * 0.2

    vmax = max(dat["Crick"].max(), dat["Watson"].max())
    vmax * 1.5

    chroms = list(sorted(set(dat["Chrom"])))

    rows = []

    for chrom in chroms:
        d = dat[dat["Chrom"] == chrom]
        ys1 = d["Crick"].values
        ys2 = d["Watson"].values
        xs = np.arange(len(ys1))
        ps1 = get_sharp(ys1, cutoff1, cutoff1)
        ps2 = get_sharp(ys2, cutoff1, cutoff1)
        cwr = np.log2(np.divide(sum(ys1), sum(ys2)))
        rows.append([chrom, cwr, len(ps1), len(ps2)])
        plt.figure(figsize=(12, 3))
        plt.bar(xs, ys1, width=1)
        plt.bar(xs, -ys2, width=1)
        for p in ps1:
            plt.axvline(p - 0.1, color="C0", ls="--", lw=1)
        for p in ps2:
            plt.axvline(p + 0.1, color="C1", ls="--", lw=1)
        plt.xlim(min(xs) - 1, max(xs) + 1)
        plt.ylim(-vmax, vmax)
        plt.tight_layout()
        plt.savefig(outdir + "/barplot.%s.png" % chrom, dpi=300)
        plt.close()
    
    dat = pd.DataFrame(rows)
    dat.columns = ["Chrom", "CWR", "Sharp.Crick", "Sharp.Watson"]
    dat.to_csv(outdir + "/sharps.tsv", sep="\t", index=False)
    
    
if __name__ == '__main__':
    main()
    