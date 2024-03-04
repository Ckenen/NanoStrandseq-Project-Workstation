#!/usr/bin/env python
import sys
import re
import json
from collections import defaultdict
import multiprocessing as mp
import numpy as np
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import ShiftLoader

def load_snps(path, chrom=None, het_only=False, hom_only=None):
    
    snps = []
    
    with pysam.VariantFile(path) as f:
        
        sample = f.header.samples[0]
        
        try:
            
            for record in f.fetch(chrom):
                
                gt = record.samples[sample]["GT"]
                try:
                    ps = record.samples[sample]["PS"]
                except Exception:
                    ps = None
                a1 = record.alleles[gt[0]]
                a2 = record.alleles[gt[1]]
                
                if len(a1) > 1 or len(a2) > 1:
                    continue
                
                if het_only and a1 == a1:
                    continue
                
                if hom_only and a1 != a2:
                    continue
                
                snp = GRange(chrom=record.chrom, start=record.start, end=record.stop, name="SNP")
                snp.allele1 = a1
                snp.allele2 = a2
                snp.ps = ps
                
                if len(snps) == 0:
                    snps.append(snp)
                elif snp.start > snps[-1].start:
                    snps.append(snp)
                else:
                    pass
                    # print(chrom, snps[-1].start, snp.start)
                
        except Exception:
            pass
        
    snps.sort()
    
    return snps


def load_regions(path, chrom):
    regions = []
    with pysam.TabixFile(path) as f:
        try:
            for line in f.fetch(chrom):
                row = line.strip("\t").split("\t")
                chrom, start, end = row[:3]
                obj = GRange(chrom=chrom, start=int(start), end=int(end))
                regions.append(obj)
        except Exception:
            pass
    return regions


def filter_snps(snps, regions):
    loader = ShiftLoader(regions)
    snps1 = []
    for het in snps:
        if len(list(loader.fetch(obj=het))) > 0:
            snps1.append(het)
    return snps1


def load_regions(path, chrom):
    regions = []
    with pysam.TabixFile(path) as f:
        try:
            for line in f.fetch(chrom):
                row = line.strip("\t").split("\t")
                chrom, start, end = row[:3]
                obj = GRange(chrom=chrom, start=int(start), end=int(end))
                regions.append(obj)
        except Exception:
            pass
    return regions


def worker(f_reference, f_query, f_bed, chrom):
    snps1 = load_snps(f_reference, chrom)
    snps2 = load_snps(f_query, chrom)

    regions =load_regions(f_bed, chrom)
    snps1 = filter_snps(snps1, regions)
    snps2 = filter_snps(snps2, regions)
    
    total = 0
    counter = defaultdict(int)
    for snp1, snp2 in zip(snps1, snps2):
        assert snp1.chrom == snp2.chrom
        assert snp1.start == snp2.start
        k = str(snp2.ps)
        if k not in counter:
            counter[k] = [0, 0]
        if snp1.allele1 != snp1.allele2 and snp1.ps == "PATMAT":
            total += 1
            if snp1.allele1 == snp2.allele1 and snp1.allele2 == snp2.allele2:
                counter[k][0] += 1
            elif snp1.allele1 == snp2.allele2 and snp1.allele2 == snp2.allele1:
                counter[k][1] += 1
            else:
                assert False
                
    items = list(sorted(counter.items(), key=lambda item: sum(item[1]), reverse=True))
    if len(items) > 0:
        v1, v2 = items[0][1]
    else:
        v1, v2 = 0, 0

    d = dict()
    d["Chrom"] = chrom
    d["HET_SNP_Reference"] = total
    d["HET_SNP_Query"] = v1 + v2
    d["HET_SNP_Overlap"] = max(v1, v2)
    d["Phasing_Recall"] = np.divide(d["HET_SNP_Overlap"], d["HET_SNP_Reference"])
    d["Phasing_Precision"] = np.divide(d["HET_SNP_Overlap"], d["HET_SNP_Query"])
    print(d)
    return d


if __name__ == "__main__":
    f_ref, f_que, f_bed, threads, outfile = sys.argv[1:]
    threads = int(threads)
    
    chroms = []
    with pysam.VariantFile(f_ref) as f:
        chroms.extend(f.header.contigs)
    with pysam.VariantFile(f_que) as f:
        chroms.extend(f.header.contigs)
    chroms = list(sorted(set(chroms)))
        
    results = []
    pool = mp.Pool(threads)
    for chrom in chroms:
        r = pool.apply_async(worker, (f_ref, f_que, f_bed, chrom))
        results.append(r)
    pool.close()
    pool.join()
    
    all_results = [r.get() for r in results]
    
    results = list(filter(lambda item: re.match("^chr[0-9]+$", item["Chrom"]), all_results))
    
    d = dict()
    d["HET_SNP_Reference"] = sum([r["HET_SNP_Reference"] for r in results])
    d["HET_SNP_Query"] = sum([r["HET_SNP_Query"] for r in results])
    d["HET_SNP_Overlap"] = sum([r["HET_SNP_Overlap"] for r in results])
    d["Phasing_Recall"] = np.divide(d["HET_SNP_Overlap"], d["HET_SNP_Reference"])
    d["Phasing_Precision"] = np.divide(d["HET_SNP_Overlap"], d["HET_SNP_Query"])
    d["Chromosomes"] = all_results
    
    with open(outfile, "w+") as fw:
        json.dump(d, fw, indent=4)
    
    