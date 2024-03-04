#!/usr/bin/env python
import sys
import pysam
from pyBioInfo.IO.File import BedFile
from pyBioInfo.Utils import ShiftLoader


def main():
    infile, blacklist, outfile = sys.argv[1:]

    with BedFile(blacklist) as f:
        regions = [x for x in f]

    loader = ShiftLoader(regions)
    with pysam.AlignmentFile(infile) as f, pysam.AlignmentFile(outfile, "wb", f) as fw:
        for chrom in sorted(f.references):
            for segment in f.fetch(chrom):
                keep = True
                start = segment.reference_start
                end = segment.reference_end
                for region in loader.fetch(chrom=chrom, start=start, end=end):
                    keep = False
                    break
                if keep:
                    fw.write(segment)


if __name__ == '__main__':
    main()
