#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    infiles = sys.argv[1:-1]
    outfile = sys.argv[-1]
    
    run = outfile.split("/")[-1][:-4]
    
    markers = ["o", "s", "P", "X", "*", "p", "D", "<", ">", "v"] * 2
    fig, axs = plt.subplots(1, 2, figsize=(20, max(len(infiles) / 12, 10)))
    plt.sca(axs[0])
    plt.title(run)
    i = 0
    for path in infiles:
        cell = path.split("/")[-1][:-4]
        d = pd.read_csv(path, sep="\t")
        # d = d[d["Sample"] > 0]
        xs = d["Sample"].values
        xs = np.log10(xs)
        ys = 1 - d["UniqRatio"].values
        ys = ys * 100
        plt.plot(xs, ys, label=cell, lw=1, marker=markers[int(i/20)], color="C%d" % (i%10))
        i += 1
    plt.ylim(0, 100)
    plt.xlabel("Selected alignments (log$_{10}$)")
    plt.ylabel("PCR duplication (%)")
    plt.grid()
    plt.legend(markerscale=1, ncol=3, bbox_to_anchor=(1, 1), loc="upper left")
    axs[1].set_visible(False)
    plt.savefig(outfile, dpi=300)
    
    
if __name__ == '__main__':
    main()
    