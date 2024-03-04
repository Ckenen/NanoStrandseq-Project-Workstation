#!/usr/bin/env python
import sys
import re
import pandas as pd
import pysam


def parse_cigar_string(s):
    array = []
    while len(s) > 0:
        ret = re.search("^[0-9]+[MIDSH]", s)
        i1, i2, = ret.span()
        s1 = s[:i2]
        k, v = s1[-1], int(s1[:-1])
        s = s[i2:]
        array.append([k ,v])
    return array


def get_mapped_region(cigars, start):
    end = start
    for cigar in cigars:
        if cigar[0] in ["M", "D", "N"]:
            end += cigar[1]
    return start, end


def main():
    f_bam, f_out = sys.argv[1:]
    
    rows = []
    with pysam.AlignmentFile(f_bam) as f:
        for s in f:
            if s.is_secondary:
                continue
            if s.is_supplementary:
                continue
            if s.mapping_quality < 30:
                continue
            chrom = s.reference_name
            chrom_length = f.get_reference_length(chrom)
            start = s.reference_start
            end = s.reference_end
            name = s.query_name
            contig_length = s.query_length
            strand = "+" if s.is_forward else "-"
            row = [chrom, chrom_length, start, end, end - start, name, contig_length, s.mapping_quality, strand]
            rows.append(row)
            
            if s.has_tag("SA"):
                for sa in s.get_tag("SA").split(";"):
                    if sa == "":
                        continue
                    sa = sa.split(",")
                    chrom1, start1, strand1, cigars1, mapq1, nm1 = sa
                    chrom_length1 = f.get_reference_length(chrom1)
                    cigars1 = parse_cigar_string(cigars1)
                    start1, end1 = get_mapped_region(cigars1, int(start1))
                    row = [chrom1, chrom_length1, start1, end1, end1 - start1, name, contig_length, mapq1, strand1]
                    if chrom1 == chrom:
                        rows.append(row)
    mapped_info = pd.DataFrame(rows)
    mapped_info.columns = ["Chrom", "ChromLength", "MappedStart", "MappedEnd", "MappedLength", 
                           "Contig", "ContigLength", "MappingQuality", "Strand"]
    mapped_info.to_csv(f_out, sep="\t", index=False)
        
    
if __name__ == "__main__":
    main()