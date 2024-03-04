#!/usr/bin/env python
import sys
from collections import defaultdict
import re
import numpy as np
import pandas as pd
import pysam
from pyBioInfo.IO.File import BedFile


def main():
    # input.bam Ns.bed blackhole.bed output.tsv
    infile1, infile2, infile3, outfile = sys.argv[1:]
    
    bin_width = 1000000
    
    blanks = defaultdict(list)
    with BedFile(infile2) as f:
        for obj in f:
            blanks[obj.chrom].append([obj.start, obj.end])
            
    blackholes = defaultdict(list)
    with BedFile(infile3) as f:
        for obj in f:
            blackholes[obj.chrom].append([obj.start, obj.end])
            
    d1 = defaultdict(int)
    for chrom, items in blanks.items():
        for start, end in items:
            i1 = int(start / bin_width)
            i2 = int((end - 1) / bin_width) + 1
            for i in range(i1, i2):
                start1, end1 = i * bin_width, (i + 1) * bin_width
                start2, end2 = max(start, start1), min(end, end1)
                d1[(chrom, i)] += (end2 - start2)
                
    d2 = defaultdict(int)
    for chrom, items in blackholes.items():
        for start, end in items:
            i1 = int(start / bin_width)
            i2 = int((end - 1) / bin_width) + 1
            for i in range(i1, i2):
                start1, end1 = i * bin_width, (i + 1) * bin_width
                start2, end2 = max(start, start1), min(end, end1)
                d2[(chrom, i)] += (end2 - start2)
            
    chroms = []
    lengths = dict() 
    d3 = defaultdict(int)
    with pysam.AlignmentFile(infile1) as f:
        for chrom in f.references:
            length = f.get_reference_length(chrom)
            if re.match("^chr([0-9]+|[X])$", chrom) is None:
                continue
            chroms.append(chrom)
            lengths[chrom] = length            
        for s in f:
            if s.is_duplicate:
                continue
            if s.get_tag("BH") == "Y":
                continue
            chrom = s.reference_name
            bi = int(s.reference_start / bin_width)
            cw = "-" if s.is_reverse else "+"
            d3[(chrom, bi, cw)] += 1
            
    rows = []
    for chrom in chroms:
        length = lengths[chrom]
        nbin = int(length / bin_width)
        if length % bin_width > 0:
            nbin += 1
        for bi in np.arange(nbin):
            start = bi * bin_width
            end = min((bi + 1) * bin_width, length)
            c = d3[(chrom, bi, "+")]
            w = d3[(chrom, bi, "-")]
            blank = d1[(chrom, bi)]
            blackhole = d2[(chrom, bi)]
            rows.append([chrom, length, bi, start, end, end - start, c, w, blank, blackhole])
            
    dat = pd.DataFrame(rows)
    dat.columns = ["Chrom", "Length", "Bin", "Start", "End", "Width", "Crick", "Watson", "Blank", "Blackhole"]
    dat["Base"] = dat["Width"] - dat["Blank"] - dat["Blackhole"]
    
    dat.to_csv(outfile, sep="\t", index=False)    
    
    
if __name__ == '__main__':
    main()
    