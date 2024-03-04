#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    infile, outfile = sys.argv[1:]
    
    dat = pd.read_csv(infile, sep="\t")
    
    tmp = dat[dat["Autosomal"]]
    rpm_mean = tmp["RPM"].mean()
    rpm_std = tmp["RPM"].std()
    std_ratio = rpm_std / rpm_mean

    xs = np.arange(len(dat))
    ys = dat["RPM"]
    xticks = dat["Chrom"]
    colors = ["C0" if x else "C1" for x in dat["Autosomal"]]
    
    plt.figure(figsize=(8, 3))
    plt.title("File: %s, Mean: %.2f, Std: %.2f, Ratio: %.4f" % (outfile.split("/")[-1], rpm_mean, rpm_std, std_ratio))
    plt.bar(xs, ys, edgecolor="black", color=colors)
    plt.axhline(rpm_mean, color="red")
    plt.xlim(min(xs) - 0.5, max(xs) + 0.5)
    plt.xticks(xs, xticks, rotation=45)
    plt.ylabel("Reads / Million Base Range")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

    
if __name__ == '__main__':
    main()
    
    
    