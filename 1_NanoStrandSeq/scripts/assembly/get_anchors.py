#!/usr/bin/env python
import sys
import gzip


def main():
    infile = sys.argv[1]
    with gzip.open(infile, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            row = line.strip("\n").split("\t")
            fmts = row[8].split(":")
            try:
                chrom = row[0]
                pos = int(row[1])
                ref = row[3]
                alts = row[4].split(",")
                if len(ref) != 1:
                    continue
                if len(alts) != 1:
                    continue
                alt = alts[0]
                if len(alt) != 1:
                    continue
                idx = fmts.index("GT")
                count = 0
                for sample in row[9:]:
                    sample = sample.split(":")
                    gt = sample[idx]
                    if gt == "0/1" or gt == "1/0" or gt == "0|1" or gt == "1|0":
                        count += 1
                if count == 3:
                    print(chrom, pos - 1, pos, "%s/%s" % (ref, alt), sep="\t")
            except ValueError:
                continue
            
            
if __name__ == '__main__':
    main()
    
