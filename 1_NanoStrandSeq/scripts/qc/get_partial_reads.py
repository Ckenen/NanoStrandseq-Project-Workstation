#!/usr/bin/env python
import sys
import gzip
from Bio import SeqIO


def main():
    infile, outfile = sys.argv[1:]
    
    w = 200
    w2 = w * 2
    n = 0
    with gzip.open(infile, "rt") as f, open(outfile, "w+") as fw:
        for r in SeqIO.parse(f, "fastq"):
            if len(r) >= w2:
                SeqIO.write(r[:w] + r[-w:], fw, "fastq")
                n += 1
                if n >= 100000:
                    break

                
if __name__ == "__main__":
    main()