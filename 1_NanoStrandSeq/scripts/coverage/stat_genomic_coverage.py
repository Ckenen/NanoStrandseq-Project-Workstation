#!/usr/bin/env python
import sys
import os
import random
import pandas as pd
import pysam


def load_regions(path):
    lengths = dict()
    data = dict()
    with pysam.AlignmentFile(path) as sam:
        for chrom, length in zip(sam.references, sam.lengths):
            lengths[chrom] = length
        for segment in sam:
            if segment.is_unmapped:
                continue
            chrom = segment.reference_name
            start = segment.reference_start
            end = segment.reference_end
            if chrom not in data:
                data[chrom] = list()
            data[chrom].append([start, end])
    return data, lengths


def merge_regions(regions):
    regions = regions.copy()
    i = 1
    while i < len(regions):
        assert regions[i][0] >= regions[i - 1][0]
        if regions[i][0] <= regions[i - 1][1]:
            regions[i - 1][1] = max(regions[i - 1][1], regions[i][1])
            regions.pop(i)
        else:
            i += 1
    return regions


def get_coverage(percentage, n, data):
    random.seed(percentage + n)
    coverages = dict()
    total = 0  # total regions
    sampling = 0  # sampling regions
    for chrom in data:
        regions = []
        for region in data[chrom]:
            total += 1
            if random.random() < percentage:
                regions.append(region)
                sampling += 1
        regions = merge_regions(regions)
        coverages[chrom] = sum([end - start for start, end in regions])
    return sampling, total, coverages


def main():
    infile, outdir = sys.argv[1:]

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    data, lengths = load_regions(infile)
    percentages = [0.01, 0.05,
                   0.1, 0.2, 0.3, 0.4, 0.5,
                   0.6, 0.7, 0.8, 0.9, 1.0]
    # percentages = [0.01, 0.05, 0.1]
    # repeat = 2
    repeat = 1
    
    results = []
    for p in percentages:
        for n in range(repeat):
            print("Try to sampling %.2f%% reads, %d times." % (p * 100, n))
            sampling, total, coverages = get_coverage(p, n, data)
            print("Exactly sampling %d (%.2f%%) reads from %d reads." %
                  (sampling, sampling * 100 / total, total))
            rows = []
            for chrom, length in lengths.items():
                cov = coverages.get(chrom, 0)
                rows.append([chrom, cov, length])
            dat = pd.DataFrame(rows)
            dat.columns = ["chrom", "coverage", "length"]
            dat["percentage"] = dat["coverage"] / dat["length"]
            genomic_coverage = dat["coverage"].sum() / dat["length"].sum()
            print("Genomic coverage is %.2f%%" % (genomic_coverage * 100))
            results.append([p, n, sampling, total, genomic_coverage, dat])
            if p == 1:
                break
    
    with open(outdir + "/genomic_coverage.tsv", "w+") as fw:
        fw.write("percentage\trepeat\tsampling\ttotal\tgenomic_coverage\n")
        for result in results:
            line = "\t".join(map(str, result[:5]))
            fw.write(line + "\n")
    
    for result in results:
        path = outdir + "/sampling_%s_%s.tsv" % (result[0], result[1])
        dat.to_csv(path, sep="\t", index=False)
    
    print("Finished!")


if __name__ == '__main__':
    main()
