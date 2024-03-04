#!/usr/bin/env python
import sys
from collections import defaultdict
import pandas as pd
import re
from pyBioInfo.IO.File import BedFile, BamFile


def main():
    infile1, infile2, outfile = sys.argv[1:]
    
    counter = defaultdict(int)
    chrom_lengths = dict()
    with BamFile(infile1) as f:
        sam_header = f.handle.header.as_dict()
        for item in sam_header["SQ"]:
            chrom_lengths[item["SN"]] = item["LN"]
        for align in f:
            chrom = align.chrom
            counter[chrom] += 1
            
    n_regions = defaultdict(list)
    with BedFile(infile2) as f:
        for region in f:
            n_regions[region.chrom].append([region.start, region.end])
            
    total_reads = sum(counter.values())
    rows = []
    for chrom, length in chrom_lengths.items():
        ns = sum([y - x for x, y in n_regions[chrom]])
        nsr = ns / length
        count = counter[chrom]
        cr = count / total_reads
        rows.append([chrom, length, ns, nsr, count, cr])
    dat = pd.DataFrame(rows)
    dat.columns = ["Chrom", "Length", "Ns", "NRatio", "Reads", "ReadRatio"]
    dat["RPM"] = dat["Reads"] * 1e6 / (dat["Length"] - dat["Ns"])
    dat["Autosomal"] = [re.match("^chr[0-9]+$", chrom) is not None for chrom in dat["Chrom"]]
    dat.to_csv(outfile, sep="\t", index=False)
    
    
if __name__ == '__main__':
    main()
    
    
    