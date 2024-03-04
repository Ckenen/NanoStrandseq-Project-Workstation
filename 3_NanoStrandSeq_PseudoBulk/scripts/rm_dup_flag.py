#!/usr/bin/env python
import sys
import pysam


def main():
    infile, outfile = sys.argv[1:]

    with pysam.AlignmentFile(infile) as f, pysam.AlignmentFile(outfile, "wb", f) as fw:
        for segment in f.fetch(until_eof=True):
            segment.flag = segment.flag & 0b101111111111
            fw.write(segment)


if __name__ == "__main__":
    main()