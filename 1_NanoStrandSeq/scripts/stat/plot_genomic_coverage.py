#!/usr/bin/env python
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    infiles = sys.argv[1:-1]
    outfile = sys.argv[-1]
    
    run = outfile.split("/")[-1][:-4]
    ncol = 3
    nrow = len(infiles) / ncol
    if len(infiles) % ncol > 0:
        nrow += 1
    fig, axs = plt.subplots(1, 2, figsize=(20, max(nrow / 4, 10)))
    plt.sca(axs[0])
    plt.title(run)
    markers = ["o", "s", "P", "X", "*", "p", "D", "<", ">", "v"] * 2
    i = 0
    vs = []
    for path in infiles:
        cell = path.split("/")[-1][:-4]
        if not os.path.exists(path):
            continue
        d = pd.read_csv(path, sep="\t")
        # if d["TotalReads"].values[0] < 200000:
        #     continue
        xs = d["Sample"].values / 1e6
        ys = d["Coverage"].values * 100
        marker = markers[int(i/20)]
        color = "C%d" % (i%10)
        plt.plot(xs, ys, lw=0.5, color=color, marker=marker, label=cell, markersize=3)
        i += 1
        vs.append(d[d["Percentage"] == 1]["Coverage"].values[0])
    mean = np.mean(vs) * 100
    median = np.median(vs) * 100
    x1, x2 = plt.gca().get_xlim()
    y1, y2 = plt.gca().get_ylim()
    w = x2 - x1
    h = y2 - y1
    plt.text(x1 + w * 0.05, y1 + h * 0.8, "mean=%.2f%%" % mean)
    plt.text(x1 + w * 0.05, y1 + h * 0.75, "median=%.2f%%" % median)
    plt.xlabel("Selected reads (1e6)")
    plt.ylabel("Genomic coverage (%)")
    plt.grid()
    plt.legend(markerscale=1, ncol=3, bbox_to_anchor=(1, 1), loc="upper left")
    axs[1].set_visible(False)
    plt.savefig(outfile, dpi=300)


if __name__ == '__main__':
    main()
    
