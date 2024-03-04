#!/usr/bin/env python
import sys
import os
import re
import numpy as np
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams["font.family"] = "arial"
import matplotlib.pyplot as plt
import pysam


def retain_cc_or_ww(flags):
    tmp1 = []
    i1 = None
    for i2, v in enumerate(flags):
        if v == 1:
            if i1 is None:
                i1 = i2
        else:
            if i1 is not None:
                tmp1.append([i1, i2])
                i1 = None
    if i1 is not None:
        tmp1.append([i1, len(flags)])

    tmp2 = []
    for i1, i2 in tmp1:
        if i2 - i1 >= 4:
            tmp2.append([i1, i2])

    tmp3 = []
    flags1 = np.zeros(len(flags), dtype=np.int) # is WC
    for i1, i2 in tmp2:
        i1 = max(0, i1 - 3)
        i2 = min(i2 + 3, len(flags))
        for i3 in range(i1, i2):
            flags1[i3] = 1
        tmp3.append([i1, i2])

    tmp4 = [] # retain CC/WW
    i1 = None
    for i2, v in enumerate(flags1):
        if v == 0:
            if i1 is None:
                i1 = i2
        else:
            if i1 is not None:
                tmp4.append([i1, i2])
                i1 = None
    if i1 is not None:
        tmp4.append([i1, len(flags1)])
    return tmp4


def main():
    infile, outdir = sys.argv[1:]

    cell = infile.split("/")[-1][:-4]
    width = 1e6
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    plotdir = outdir + "/plots"
    if not os.path.exists(plotdir):
        os.mkdir(plotdir)

    chroms = []
    lengths = dict()
    data = dict()
    with pysam.AlignmentFile(infile) as f:
        for chrom in f.references:
            chroms.append(chrom)
            length = f.get_reference_length(chrom)
            lengths[chrom] = length
            nbin = int(length / width)
            if length % width > 0:
                nbin += 1
            counts1 = np.zeros(nbin, dtype=np.int) # crick
            counts2 = np.zeros(nbin, dtype=np.int) # watson
            for segment in f.fetch(chrom):
                if segment.is_reverse:
                    counts = counts2
                else:
                    counts = counts1
                if segment.is_duplicate:
                    continue
                if segment.get_tag("BH") == "Y":
                    continue
                counts[int(segment.reference_start / width)] += 1
            data[chrom] = [counts1, counts2]

    max_length = max(lengths.values())
    max_nbin = int(max_length / width)
    if max_length % width > 0:
        max_nbin += 1

    rows = []
    for chrom in chroms:
        if re.match("^chr([0-9]+|[XY])$", chrom) is None:
            continue

        length = lengths[chrom]
        counts1, counts2 = data[chrom]
        xlim = max(max(max(counts1), max(counts2)) * 1.2, 10)
        ys = np.arange(len(counts1)) + 0.5
        
        fig, axs = plt.subplots(1, 3, figsize=(6, 6), sharex=True, sharey=True)
        title = "%s.%s" % (cell, chrom)
        plt.suptitle(title)
        
        ax = axs[0]
        plt.sca(ax)
        plt.title("Raw")
        plt.barh(ys, counts1, height=1, color="C0")
        plt.barh(ys, -counts2, height=1, color="C1")
        plt.plot([0, 0], [0, len(counts1)], lw=1, color="black")
        plt.xlim(-xlim, xlim)
        plt.ylim(0, max_nbin)
        
        ax = axs[1]
        plt.sca(ax)
        plt.title("Not WC")
        plt.barh(ys, counts1, height=1, color="grey")
        plt.barh(ys, -counts2, height=1, color="grey")
        plt.plot([0, 0], [0, len(counts1)], lw=1, color="black")
        xs1, xs2 = [], []
        idxs = []
        vs = counts1 + counts2
        vs = vs[vs > 0]
        cutoff = np.median(vs) * 0.5 * 0.33
        mean = np.mean(vs)
        std = np.std(vs)
        r = std / mean
        flags = np.zeros(len(counts1))
        for i, c, w in zip(np.arange(len(counts1)), counts1, counts2):
            if c >= cutoff and w >= cutoff:
                if np.abs(np.log2(np.divide(c, w))) < 1:
                    flags[i] = 1
                    c, w = 0, 0
            xs1.append(c)
            xs2.append(w)
        xs1, xs2 = np.array(xs1), np.array(xs2)
        plt.barh(ys, xs1, height=1, color="C0")
        plt.barh(ys, -xs2, height=1, color="C1")
        
        tmp4 = retain_cc_or_ww(flags)
        y1, y2 = None, None
        t = "X"
        log2cwr = 0.0
        if len(tmp4) > 0:
            tmp4 = list(sorted(tmp4, key=lambda item: item[1] - item[0]))
            item = tmp4[-1]
            if item[1] - item[0] > len(vs) * 0.33:
                y1, y2 = item
                c = np.sum(counts1[y1:y2])
                w = np.sum(counts2[y1:y2])
                log2cwr = np.log2(np.divide(c, w))
                cutoff = max(3, (item[1] - item[0] - 40) * 0.0025 + 3)
                if log2cwr > cutoff:
                    t = "CC"
                elif log2cwr < -cutoff:
                    t = "WW"
                rows.append([chrom, cell, y1, y2, t, log2cwr])
        
        ax = axs[2]
        plt.sca(ax)
        plt.title("%s [%.2f]" % (t, log2cwr))
        plt.barh(ys, counts1, height=1, color="grey")
        plt.barh(ys, -counts2, height=1, color="grey")
        if y1 is not None:
            ysA = []
            xsC = []
            xsW = []
            for y in range(y1, y2):
                ysA.append(y)
                xsC.append(counts1[y])
                xsW.append(counts2[y])
            ysA = np.array(ysA) + 0.5
            xsC = np.array(xsC)
            xsW = np.array(xsW)
            plt.barh(ysA, xsC, height=1, color="C0")
            plt.barh(ysA, -xsW, height=1, color="C1")
                    
        plt.savefig(outdir + "/plots/%s.pdf" % chrom, dpi=300)

    with open(outdir + "/ccww.tsv", "w+") as fw:
        for row in rows:
            line = "\t".join(map(str, row))
            fw.write(line + "\n")
            

if __name__ == "__main__":
    main()


