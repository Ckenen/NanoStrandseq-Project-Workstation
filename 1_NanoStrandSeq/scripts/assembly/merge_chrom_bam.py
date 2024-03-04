#!/usr/bin/env python
import os
import sys
import pysam


def main():
    infiles = sys.argv[1:-2]
    chrom = sys.argv[-2]
    outfile = sys.argv[-1]

    with pysam.AlignmentFile(infiles[0]) as f:
        header = f.header.as_dict()
        items = []
        for infile in infiles:
            name = os.path.basename(os.path.splitext(infile)[0])
            item = {'ID': name, 'LB': name, 'SM': name}
            items.append(item)
        header["RG"] = items

    with pysam.AlignmentFile(outfile, "wb", header=header) as fw:
        for infile in infiles:
            with pysam.AlignmentFile(infile) as f:
                name = os.path.basename(os.path.splitext(infile)[0])
                for segment in f.fetch(contig=chrom):
                    segment.set_tag("RG", name)
                    fw.write(segment)
                    
                    
if __name__ == '__main__':
    main()
    