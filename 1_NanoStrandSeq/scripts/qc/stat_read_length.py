#!/usr/bin/env python
import sys
from collections import defaultdict
import pysam


def main():
    infile = sys.argv[1]
    counter = defaultdict(int)
    with pysam.FastxFile(infile) as f:
        for read in f:
            counter[len(read.sequence)] += 1
    for k, v in sorted(counter.items()):
        print(k, v, sep="\t")
        
    
if __name__ == "__main__":
    main()