#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plot(title, data, full, outfile):
    chroms = ["chr%d" % i for i in range(1, 23)] + ["chrX"]
    chroms = list(filter(lambda item: item in data, chroms))
    nrow = 2
    ncol = 12 # int(len(chroms) / 2)
    if len(chroms) % 2 != 0:
        ncol += 1

    fig, axs = plt.subplots(nrow, ncol, figsize=(ncol * 2.5, nrow * 6))
    plt.suptitle(title)
    
    max_y = max([len(data[chrom][0]) for chrom in chroms]) * 1.1
    if full:
        max_x = max([max(max(data[chrom][0]), max(data[chrom][1]))
                    for chrom in chroms]) * 1.2
    else:
        values = []
        for chrom in chroms:
            values.extend(data[chrom][0])
            values.extend(data[chrom][1])
        max_x = np.mean(values) + 2 * np.std(values)

    for i in range(len(chroms)):
        irow = int(i / ncol)
        icol = i % ncol
        chrom = chroms[i]
        xs1 = data[chrom][0]
        xs2 = data[chrom][1]
        xs2 = -xs2
        ys = np.arange(len(xs1))
        plt.sca(axs[irow][icol])
        plt.xlabel(chrom)
        plt.plot([0, 0], [0, len(ys)], lw=0.1, color="black", zorder=0)
        plt.text(max_x * 0.05, max_y * 0.98, "mean=%.2f" %
                 np.mean(xs1), color="C0")
        plt.text(max_x * 0.05, max_y * 0.95, "median=%.2f" %
                 np.median(xs1), color="C0")
        plt.text(max_x * 0.05, max_y * 0.92, "std=%.2f" %
                 np.std(xs1), color="C0")
        plt.text(-max_x, max_y * 0.98, "mean=%.2f" % np.mean(xs2), color="C1")
        plt.text(-max_x, max_y * 0.95, "median=%.2f" %
                 np.median(xs2), color="C1")
        plt.text(-max_x, max_y * 0.92, "std=%.2f" % np.std(xs2), color="C1")
        plt.plot([np.mean(xs1), np.mean(xs1)], [0, len(ys)],
                 color="black", ls="--", lw=0.1, zorder=0)
        plt.plot([np.mean(xs2), np.mean(xs2)], [0, len(ys)],
                 color="black", ls="--", lw=0.1, zorder=0)
        plt.barh(ys, xs1, height=1, color="C0")
        plt.barh(ys, xs2, height=1, color="C1")
        array = data[chrom]
        if len(array) == 6:
            xs3 = array[2]
            xs4 = array[3]
            xs5 = -array[4]
            xs6 = -array[5]
            plt.barh(ys, xs3, color="blue")
            plt.barh(ys, xs4, left=xs3, color="red")
            plt.barh(ys, xs5, color="blue")
            plt.barh(ys, xs6, left=xs5, color="red")

    # 隐藏边框
    for ax in axs.flatten():
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # ax.set_xticks([-max_x, 0, max_x])
        ax.set_yticks([])
        ax.set_xlim(-max_x, max_x)
        ax.set_ylim(-1, max_y + 1)

    # 屏蔽掉多余的axes
    n = 0
    for i in range(nrow):
        for j in range(ncol):
            n += 1
            if n > len(chroms):
                axs[i][j].set_visible(False)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)


def main():
    infile, prefix = sys.argv[1:]

    with open(infile) as f:
        data = dict()
        for line in f:
            row = line.strip("\n").split("\t")
            chrom = row[0]
            wbin = int(row[1])
            array = [np.array(list(map(float, values.split(",")))) for values in row[2:]]
            # chrom, wbin, values1, values2 = line.strip("\n").split("\t")
            # values1 = np.array(list(map(float, values1.split(","))))
            # values2 = np.array(list(map(float, values2.split(","))))
            # values1 = array[0]
            # values2 = array[1]
            # data[chrom] = [values1, values2]
            data[chrom] = array
    title = infile
    plot(title, data, True, prefix + ".full.pdf")
    plot(title, data, False, prefix + ".trimmed.pdf")


if __name__ == '__main__':
    main()
