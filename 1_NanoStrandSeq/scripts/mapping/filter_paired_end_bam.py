#!/usr/bin/env python
import sys
import re
from collections import defaultdict
import pysam

def cal_perc(v1, v2):
    if v2 == 0:
        return 0
    else:
        return v1 * 100 / v2

def main():
    infile, outfile = sys.argv[1:]

    n_total = 0
    n_filtered = 0
    n_paired = 0
    
    with pysam.AlignmentFile(infile) as f, pysam.AlignmentFile(outfile, "wb", f) as fw:
        for chrom in f.references:
            
            if re.match("^chr([0-9]+|[XY])$", chrom) is None:
                continue
            
            segments = defaultdict(list)
            for s in f.fetch(chrom):
                n_total += 1
                if s.is_secondary:
                    continue
                if s.is_supplementary:
                    continue
                if not s.is_proper_pair:
                    continue
                if s.mapping_quality < 20:
                    continue
                segments[s.query_name].append(s)
                n_filtered += 1
                
            paired_segments = []
            for k, v in segments.items():
                if len(v) != 2:
                    continue
                r1, r2 = v
                if r1.is_read2:
                    r1, r2 = r2, r1
                assert r1.is_read1 and r2.is_read2
                paired_segments.append(r1)
                paired_segments.append(r2)
                n_paired += 2
                
            for s in sorted(paired_segments, key=lambda s: s.reference_start):
                fw.write(s)
                
    print("Input reads: %d (%.2f%%)" % (n_total, cal_perc(n_total, n_total)))
    print("Output reads: %d (%.2f%%)" % (n_paired, cal_perc(n_paired, n_total)))
    
    
if __name__ == "__main__":
    main()

