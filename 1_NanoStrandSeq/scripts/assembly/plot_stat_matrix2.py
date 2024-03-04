#!/usr/bin/env python
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("AGG")
matplotlib.rcParams["font.family"] = "arial"
import matplotlib.pyplot as plt


def main():
    infile, outdir = sys.argv[1:]
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    dat = pd.read_csv(infile, sep="\t")
    ds = [dat, dat[dat["HC"] == 0], dat[dat["HC"] == 1]]
    
    for i, d in enumerate(ds):    
        for hp in [1, 2]:
            prefix = outdir + "/%s" % os.path.basename(outdir)
            if i == 0:
                prefix += ".all"
            elif i == 1:
                prefix += ".HC0"
            else:
                prefix += ".HC1"
            prefix += ".hp%d" % hp
            
            cases = np.arange(1, 7)
            confs = [True, False]
            rows = []
            for case in cases:
                row = []
                for conf in confs:
                    v = d[(d["Case%d" % hp] == case) & (d["Conf%d" % hp] == conf)]["Count"].sum()
                    row.append(v)
                rows.append(row)
            df = pd.DataFrame(rows)
            df.index = cases
            df.index.name = "Case"
            df.columns = ["Consistent", "Inconsistent"]
            df["Total"] = df["Consistent"] + df["Inconsistent"]
            total = df["Total"].sum()
            df["Total%"] = df["Total"] / total
            df["Consistent.Perc"] = df["Consistent"] / total
            df["Inconsistent.Perc"] = df["Inconsistent"] / total
            df.to_csv(prefix + ".tsv", sep="\t")
            
            xs = cases
            ys1 = df["Consistent.Perc"] * 100
            ys2 = df["Inconsistent.Perc"] * 100
            
            plt.figure(figsize=(4, 3))
            plt.title(os.path.basename(prefix))
            plt.bar(xs, ys1, edgecolor="black", width=0.6, label="Consistent")
            plt.bar(xs, ys2, edgecolor="black", width=0.6, bottom=ys1, label="Inconsistent")
            plt.xlabel("Case")
            plt.ylabel("Percentage (%)")
            plt.xlim(min(xs) - 0.5, max(xs) + 0.5)
            plt.ylim(0, 100)
            plt.legend()
            plt.tight_layout()
            plt.savefig(prefix + ".pdf", dpi=300)
            
            
if __name__ == "__main__":
    main()