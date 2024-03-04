#!/usr/bin/env python
import sys
from collections import defaultdict
import multiprocessing
import re
import numpy as np
import pyBigWig

def stat_chrom_base_depth(infile, chrom):
    step = 1000000
    counter = defaultdict(int)
    with pyBigWig.open(infile) as f:
        length = f.chroms()[chrom]
        for start in range(0, length, step):
            end = min(start + step, length)
            for v in f.values(chrom, start, end):
                if np.isnan(v):
                    v = 0
                else:
                    v = int(v)
                counter[v] += 1
    return counter

def main():
    # input.bw threads output.tsv
    infile, threads, outfile = sys.argv[1:]
    threads = int(threads)

    results = []
    pool = None
    if threads > 1:
        pool = multiprocessing.Pool(threads)
    with pyBigWig.open(infile) as f:
        for chrom in f.chroms():
            if re.match("^chr[0-9]+$", chrom) is None:
                continue
            if pool:
                results.append(pool.apply_async(stat_chrom_base_depth, (infile, chrom)))
            else:
                results.append(stat_chrom_base_depth(infile, chrom))
    if pool:
        pool.close()
        pool.join()
        results = [r.get() for r in results]
            
    counter = defaultdict(int)        
    for r in results:
        for k, v in r.items():
            counter[k] += v
    with open(outfile, "w+") as fw:
        for k, v in sorted(counter.items()):
            line = "%d\t%d" % (k, v)
            fw.write(line + "\n")


if __name__ == "__main__":
    main()