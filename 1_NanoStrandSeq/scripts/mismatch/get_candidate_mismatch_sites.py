#!/usr/bin/env python
import sys
import gzip

infile, outfile = sys.argv[1:]
with gzip.open(infile, "rt") as f, gzip.open(outfile, "wt") as fw:
    # INPUT: chrom position strand events
    for line in f:
        line = line.strip("\n")
        chrom, position, strand, events = line.split("\t")
        num = 0
        for event in events.split(";"):
            if event == "":
                continue
            ref, alt, qua = event.split(",")
            qua = int(qua)
            if qua >= 8:
                num += 1

        if num >= 2:
            # OUTPUT: chrom position strand
            line = "\t".join([chrom, position])
            fw.write(line)
            fw.write("\n")
