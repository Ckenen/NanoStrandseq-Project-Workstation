#!/usr/bin/env python
import os
import sys
import pysam


# def make_header(infiles):
#     headers = []
#     for infile in infiles:
#         with pysam.AlignmentFile(infile) as f:
#             headers.append(f.header.as_dict())
#     header = headers[0].copy()
#     for h in headers[1:]:
#         for item in h["RG"]:
#             header["RG"].append(item)
#     return header


def main():
    infiles = sys.argv[1:-2]
    chrom = sys.argv[-2]
    outdir = sys.argv[-1]

    os.mkdir(outdir)
    for infile in infiles:
        assert infile.endswith(".bam")
        outfile = "%s/%s" % (outdir, os.path.basename(infile))
        with pysam.AlignmentFile(infile) as f, pysam.AlignmentFile(outfile, "wb", f) as fw:
            for segment in f.fetch(chrom):
                fw.write(segment)


if __name__ == '__main__':
    main()
