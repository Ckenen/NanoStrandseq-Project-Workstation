#!/usr/bin/env python
import sys
import os
import gzip
import multiprocessing
import pandas as pd
import pysam


def load_ref_snvs(vcfpath, chrom, start, end):
    ref_snvs = dict()
    with pysam.AlignmentFile(vcfpath) as f:
        sample = list(f.header.samples)[0]
        for record in f.fetch(chrom, start, end):
            gt = record.samples[sample]["GT"]
            a1, a2 = record.alleles[gt[0]], record.alleles[gt[1]]
            if len(a1) != 1 or len(a2) != 1:
                continue
            ref_snvs[record.start] = [a1, a2]
    return ref_snvs


def process_bin_matrix(mpath, vcfpath, chrom, start, end):
    ref_snvs = load_ref_snvs(vcfpath, chrom, start, end)
    
    with gzip.open(mpath, "rt") as f:
        for line in f:
            pass
    

def main():
    indir, vcfpath, threads, outfile = sys.arga[1:]
    threads = int(threads)
    meta = pd.read_csv(os.path.join(indir, "meta.tsv"), sep="\t", header=None)
    pool = multiprocessing.Pool(threads)
    results = []
    for chrom, start, end, fname in meta.values:
        mpath = os.path.join(indir, fname) # matrix path
        args = (mpath, vcfpath, chrom, start, end)
        results.append(pool.apply_async(None, args))
    pool.close()
    pool.join()
    results = [r.get() for r in results]
    
    