#!/usr/bin/env python
import sys
import multiprocessing
import pandas as pd
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import BundleBuilder

def load_segments(bamfile, chrom):
    with pysam.AlignmentFile(bamfile) as f:
        for segment in f.fetch(chrom):
            obj = GRange(chrom=segment.reference_name,
                start=segment.reference_start, 
                end=segment.reference_end)
            yield obj

def execute_chromosome(bamfile, chrom):
    rows = []
    for b in BundleBuilder(load_segments(bamfile, chrom), min_capacity=10000):
        row = [b.chrom, b.start_min, b.start_max, b.end_min, b.end_max, b.count]
        rows.append(row)
    return rows

def main():
    bamfile, threads, outfile = sys.argv[1:]
    threads = int(threads)

    pool = None
    results = []
    if threads > 1:
        pool = multiprocessing.Pool(threads)
    with pysam.AlignmentFile(bamfile) as f:
        for chrom in f.references:
            if pool is None:
                results.append(execute_chromosome(bamfile, chrom))
            else:
                results.append(pool.apply_async(execute_chromosome, (bamfile, chrom)))
    if pool is not None:
        pool.close()
        pool.join()
        results = [r.get() for r in results]

    rows = []
    for rows1 in results:
        rows.extend(rows1)
    rows.sort()

    dat = pd.DataFrame(rows)
    dat.columns = ["Chrom", "StartMin", "StartMax", "EndMin", "EndMax", "Count"]
    dat.to_csv(outfile, sep="\t", index=False)
    
    # with BamFile(infile) as f:
    #     for i, bundle in enumerate(BundleBuilder(f, min_capacity=10000, keep=True)):
    #         print(bundle)
    
    
if __name__ == "__main__":
    main()