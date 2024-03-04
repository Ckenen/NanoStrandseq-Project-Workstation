#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import pysam

def main():
    infile, datadir, chrom, outfile = sys.argv[1:]
    
    bin_width = 1e6    
    
    dat = pd.read_csv(infile, sep="\t", header=None)
    dat.columns = ["Chrom", "Cell", "Index1", "Index2", "Type", "Log2CWR"]
    dat = dat[dat["Chrom"] == chrom]
    
    regions = []        
        
    for cell, idx1, idx2, log2cwr in dat[["Cell", "Index1", "Index2", "Log2CWR"]].values:
        tran = False
        cutoff = 3
        cutoff = max(3, (idx2 - idx1 - 40) * 0.0025 + 3)
        if log2cwr > cutoff:
            tran = False
        elif log2cwr < -cutoff:
            tran = True
        else:
            continue 
        run = cell.split(".")[0]
        path = datadir + "/blackhole/mark_blackhole/%s.bam" % cell
        # print(path)
        with pysam.AlignmentFile(path) as f:
            length = f.get_reference_length(chrom)
            start = idx1 * bin_width
            end = min(idx2 * bin_width, length)
            # print(start, end, length, sep="\t")
            for segment in f.fetch(chrom, start, end):
                if segment.is_duplicate:
                    continue
                if segment.get_tag("BH") == "Y":
                    continue
                if start <= segment.reference_start < end:
                    strand = "-" if segment.is_reverse else "+"
                    if tran:
                        strand = "+" if strand == "-" else "-"
                    regions.append([segment.reference_start, segment.reference_end, strand])  
    
    regions.sort()  
    with open(outfile, "w+") as fw:
        for start, end, strand in regions:
            color = "107,137,138"
            name = "C"
            if strand == "-":
                color = "248,173,97"
                name = "W"
            line = "\t".join(map(str, [chrom, start, end, name, ".", strand, start, end, color]))
            fw.write(line + "\n")
    
        
if __name__ == "__main__":
    main()