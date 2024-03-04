#!/usr/bin/env python
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam

def load_coverage(path):
    data = dict()
    wbin = 1000000
    with pysam.AlignmentFile(path) as f:
        for item in f.header.as_dict()["SQ"]:
            sn = item["SN"]
            ln = item["LN"]
            nbin = int(ln / wbin)
            if ln % wbin > 0:
                nbin += 1
            data[sn] = [np.zeros(nbin), np.zeros(nbin)]
        for segment in f:
            if segment.get_tag("XT") == "F":
                if segment.is_reverse:
                    s = 1
                else:
                    s = 0
            else:
                if segment.is_reverse:
                    s = 0
                else:
                    s = 1
            # s = 1 if segment.is_reverse else 0
            counts = data[segment.reference_name][s]
            a = int(segment.reference_start / wbin)
            b = int((segment.reference_end - 1) / wbin)
            c = b - a + 1
            w = 1 / c
            for x in range(a, b + 1):
                counts[x] += w
    return data


def main():
    infile, outfile = sys.argv[1:]
    
    data = load_coverage(infile)
                
    max_l = max([len(data[chrom][0]) for chrom in data])
    chroms = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
    # print(chroms)

    fig, axs = plt.subplots(2, 12, figsize=(12 * 1.5, 8), sharex=True, sharey=True, 
                            gridspec_kw={"wspace": 0})

    array = []
    for chrom in chroms:
        values1, values2 = data[chrom]
        array.append(values1)
        array.append(values2)
    values = np.concatenate(array, axis=0)
    width = np.mean(values) + 2 * np.std(values)
    
    for i, chrom in enumerate(chroms):
        values1, values2 = data[chrom]
        xs = np.arange(len(values1))
        ax = axs[int(i/12)][i%12]
        plt.sca(ax)
        plt.barh(xs, values1, height=1, zorder=1)
        plt.barh(xs, -values2, height=1, zorder=2)
        plt.xlim(-width, width)
        plt.ylim(-0.5, max_l + 0.5)
        plt.fill_betweenx([0, len(values1)], -width * 0.02, [width * 0.02, width * 0.02], color="lightgrey", zorder=10)
        plt.yticks([])
        plt.xticks([])
        plt.xlabel(chrom)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
    axs[-1][-1].set_visible(False)
    plt.savefig(outfile, dpi=300)
    
if __name__ == "__main__":
    main()

