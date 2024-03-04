#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd


def main():
    # bin_reads.bed blacklists.bed
    infile, outfile = sys.argv[1:]

    dat = pd.read_csv(infile, sep="\t", header=None)
    dat.columns = ["Chrom", "Start", "End", "Name", "Score", "Strand"]
    scores = dat["Score"]
    scores = scores[scores > 0]
    ln_s = np.log(scores)
    np.mean(scores), np.mean(ln_s), np.median(ln_s), np.std(ln_s)
    threshold = np.median(ln_s) + 3 * np.std(ln_s[ln_s > np.median(ln_s)])
    threshold = np.e ** threshold
    array = []
    for c, d in dat.groupby(by="Chrom"):
        flags = [False] * len(d)
        for i, v in enumerate(d["Score"] > threshold):
            if v:
                flags[i] = True
                if i > 0:
                    flags[i - 1] = True
                if i < len(d) - 1:
                    flags[i + 1] = True

        d1 = d[flags]
        array.append(d1)
    dat = pd.concat(array)
    dat.to_csv(outfile, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()