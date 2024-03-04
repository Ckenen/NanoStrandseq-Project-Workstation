#!/usr/bin/env python
import sys, glob, gzip
import multiprocessing
from collections import defaultdict
import numpy as np
import pandas as pd
from sstools.utils import BaseMatrix2


def stat_matrix2_core(path):
    counter = defaultdict(int)
    with gzip.open(path, "rt") as f:
        for line in f:
            row = BaseMatrix2.parse_line(line)
            case1 = BaseMatrix2.get_case(row[7])
            case2 = BaseMatrix2.get_case(row[8])
            conf1 = BaseMatrix2.get_confidence(row[5], row[7], case1)
            conf2 = BaseMatrix2.get_confidence(row[6], row[8], case2)
            key = (int(row[4]), case1, conf1, case2, conf2)
            counter[key] += 1
    # print(counter)
    return counter


def main():
    indir, threads, outfile = sys.argv[1:]
    threads = int(threads)
    
    results = []
    pool = multiprocessing.Pool(threads)
    for path in sorted(glob.glob(indir + "/*.matrix.gz")):
        r = pool.apply_async(stat_matrix2_core, (path,))
        results.append(r)
    pool.close()
    pool.join()
    results = [r.get() for r in results]
    
    counter = defaultdict(int)
    for r in results:
        for k, v in r.items():
            counter[k] += v
    rows = []
    for k, v in counter.items():
        row = list(k)
        row.append(v)
        rows.append(row)
    dat = pd.DataFrame(rows)
    dat.columns = ["HC", "Case1", "Conf1", "Case2", "Conf2", "Count"]
    dat.to_csv(outfile, sep="\t", index=False)
    
    
if __name__ == "__main__":
    main()
    