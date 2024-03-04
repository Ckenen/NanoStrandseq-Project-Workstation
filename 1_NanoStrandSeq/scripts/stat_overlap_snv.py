#!/usr/bin/env python
import sys
import gzip
import pyBigWig


def load(path):
    data = dict()
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            row = line.strip("\n").split("\t")
            row[1] = int(row[1])
            data[(row[0], row[1])] = row
    return data


def main():
    ref_vcf, que_vcf, bw_prefix = sys.argv[1:]

    dat1 = load(ref_vcf)
    dat2 = load(que_vcf)
    print("Reference SNVs: %d" % len(dat1))
    print("Detected SNVs: %d" % len(dat2))

    overlaps = list(set(dat1.keys() & dat2.keys()))
    print("Overlap SNVs: %d" % len(overlaps))
    type1 = 0
    type2 = 0
    genotype1 = 0
    genotype2 = 0
    for position in overlaps:
        row1 = dat1[position]
        row2 = dat2[position]
        gt1 = row1[9].split(":")[0]
        gt2 = row2[9].split(":")[0]
        if row1[3] == row2[3] and row1[4] == row2[4]:
            type1 += 1
            if gt1 == gt2:
                genotype1 += 1
            else:
                genotype2 += 1
        else:
            type2 += 1
    print("Ref/Alt concordant: %d" % type1)
    print("Genotype concordant: %d" % genotype1)
    print("Genotype discordant: %d" % genotype2)
    print("Ref/Alt discordant: %d" % type2)

    bw1 = pyBigWig.open(bw_prefix + ".+.bw")
    bw2 = pyBigWig.open(bw_prefix + ".-.bw")
    positions = list(dat1.keys() - dat2.keys())
    positions.sort()
    count1 = 0
    count2 = 0
    for chrom, end in positions:
        c1 = bw1.values(chrom, end - 1, end)[0]
        c2 = bw2.values(chrom, end - 1, end)[0]
        cov = c1 + c2
        if cov > 0:
            count1 += 1
        else:
            count2 += 1
    bw1.close()
    bw2.close()
    print("Ref only SNVs with coverage: %d" % count1)
    print("Ref only SNVs without coverage: %d" % count2)

    positions = list(dat2.keys() - dat1.keys())
    print("Que only SNVs: %d" % len(positions))


if __name__ == '__main__':
    main()
