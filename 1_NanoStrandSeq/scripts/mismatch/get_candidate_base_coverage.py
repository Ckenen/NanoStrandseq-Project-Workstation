#!/usr/bin/env python
import sys
import gzip
import pyBigWig
import numpy as np

def get_value(bw, chrom, start):
    try:
        c = bw.values(chrom, start, start + 1)[0]
        if np.isnan(c):
            c = 0
        else:
            c = int(c)
    except RuntimeError:
        c = 0
    return c

def main():
    in_tsv, in_bw1, in_bw2, out_tsv = sys.argv[1:]

    bw1 = pyBigWig.open(in_bw1)
    bw2 = pyBigWig.open(in_bw2)
    with gzip.open(in_tsv, "rt") as f, gzip.open(out_tsv, "wt") as fw:
        for line in f:
            chrom, position = line.strip("\n").split("\t")
            start = int(position)
            c1 = get_value(bw1, chrom, start)
            c2 = get_value(bw2, chrom, start)
            fw.write("%d|%d\n" % (c1, c2))
    bw1.close()
    bw2.close()


if __name__ == "__main__":
    main()
