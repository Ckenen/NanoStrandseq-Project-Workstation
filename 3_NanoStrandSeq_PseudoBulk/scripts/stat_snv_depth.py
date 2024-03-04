#!/usr/bin/env python
import sys
from collections import defaultdict
import json
import re
import multiprocessing
import numpy as np
import pysam
import pyBigWig
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import ShiftLoader


def load_benchmark_regions(path, chrom):
    regions = []
    with pysam.TabixFile(path) as f:
        for line in f.fetch(chrom):
            chrom, start, end = line.strip("\n").split("\t")
            obj = GRange(chrom=chrom, start=int(start), end=int(end))
            regions.append(obj)
    return regions


def load_snvs(path, chrom):
    snvs = []
    with pysam.VariantFile(path) as f:
        name = list(f.header.samples)[0]
        for record in f.fetch(chrom):
            gt = record.samples[name]["GT"]
            a1 = record.alleles[gt[0]]
            a2 = record.alleles[gt[1]]
            if len(a1) > 1 or len(a2) > 1:
                continue
            # if a1 == a2:
            #     continue
            end = record.pos
            start = end - 1
            obj = GRange(chrom=chrom, start=start, end=end, name="%s|%s" % (a1, a2))
            obj.a1 = a1
            obj.a2 = a2
            snvs.append(obj)
    return snvs


def filter_benchmark_snvs(snvs, regions):
    loader = ShiftLoader(regions)
    snvs1 = []
    for snv in snvs:
        if len(list(loader.fetch(obj=snv))) > 0:
            snvs1.append(snv)
    return snvs1


def get_depth(snvs, path):
    bw = pyBigWig.open(path)
    for snv in snvs:
        v = bw.values(snv.chrom, snv.start, snv.end)[0]
        if np.isnan(v):
            v = 0
        else:
            v = int(v)
        snv.score = v
    bw.close()


def stat_chrom_snv_depth(vcffile, bedfile, bwfile, chrom):
    regions = load_benchmark_regions(bedfile, chrom)
    snvs = load_snvs(vcffile, chrom)
    snvs = filter_benchmark_snvs(snvs, regions)
    get_depth(snvs, bwfile)
    counter1 = defaultdict(int)
    counter2 = defaultdict(int)
    for snv in snvs:
        if snv.a1 == snv.a2:
            counter1[snv.score] += 1
        else:
            counter2[snv.score] += 1
    return counter1, counter2


def main():
    vcffile, bedfile, bwfile, threads, outfile = sys.argv[1:]
    threads = int(threads)
    
    results = []
    pool = None
    if threads > 1:
        pool = multiprocessing.Pool(threads)
    with pysam.VariantFile(vcffile) as f:
        for chrom in f.header.contigs:
            if re.match("^chr[0-9]+$", chrom) is None:
                continue
            print(chrom)
            args = (vcffile, bedfile, bwfile, chrom)
            if pool:
                r = pool.apply_async(stat_chrom_snv_depth, args)
            else:
                r = stat_chrom_snv_depth(*args)
            results.append(r)
    if pool:
        pool.close()
        pool.join()
        results = [r.get() for r in results]
        
    counter1 = defaultdict(int)
    counter2 = defaultdict(int)
    for c1, c2 in results:
        for k, v in c1.items():
            counter1[k] += v
        for k, v in c2.items():
            counter2[k] += v
            
    data = {"HOM": counter1, "HET": counter2}
    with open(outfile, "w+") as fw:
        json.dump(data, fw)
    

if __name__ == "__main__":
    main()




