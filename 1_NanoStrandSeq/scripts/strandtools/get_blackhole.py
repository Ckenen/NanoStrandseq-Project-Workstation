#!/usr/bin/env python
import sys
import numpy as np
import pyBigWig


THRESHOLD = 10

def main():
    infile = sys.argv[1]
    
    bw = pyBigWig.open(infile)
    for chrom, length in bw.chroms().items():
        values = bw.values(chrom, 0, length)
        values = np.nan_to_num(values, 0)
        start = None
        regions = []
        for i, v in enumerate(values):
            if v >= THRESHOLD:
                if start is None:
                    start = i
            else:
                if start is not None:
                    regions.append([start, i])
                    start = None
        if start is not None:
            regions.append([start, i])

        for start, end in regions:
            print(chrom, start, end, sep="\t")
    
    bw.close()
    
if __name__ == '__main__':
    main()
    