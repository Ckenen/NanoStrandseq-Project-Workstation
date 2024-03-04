#!/usr/bin/env python
import sys
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam


def main():
    infile, outfile = sys.argv[1:]
    
    
    vs = []
    with pysam.AlignmentFile(infile) as f:
        for segment in f:
            dn = segment.get_tag("DN")
            ds = segment.get_tag("DS")
            di = segment.get_tag("DI")
            if di != 0:
                continue
            vs.append(ds)
            
    counter = Counter(vs)
    total = len(vs)
    with open(outfile, "w+") as fw:
        fw.write("ReadPerDuplicate\tNumber\tTotal\tRatio\n")
        for k in sorted(counter.keys()):
            v = counter[k]
            r = v / total
            fw.write("\t".join(map(str, [k, v, total, r])) + "\n")
    
    prefix = outfile
    if outfile.endswith(".tsv"):
        prefix = outfile[:-4]
        
    xs = np.arange(11)
    ys = np.zeros(len(xs))

    for k, v in counter.items():
        i = min(k, len(xs) - 1)
        ys[i] += v
    ys = ys * 100 / sum(ys)
    plt.figure(figsize=(4, 3))
    plt.bar(xs, ys)
    plt.xticks(xs)
    plt.xlabel("Reads per duplicate set")
    plt.ylabel("Percentage (%)")
    plt.grid(axis="y", ls="--")
    plt.tight_layout()
    plt.savefig(prefix + ".png", dpi=300)
    plt.close()
    
        
if __name__ == '__main__':
    main()
    

