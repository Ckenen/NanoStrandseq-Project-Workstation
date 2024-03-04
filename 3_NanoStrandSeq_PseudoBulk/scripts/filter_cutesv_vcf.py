#!/usr/bin/env python
import sys
import re

for line in sys.stdin:
    if line.startswith("#"):
        if line.startswith("##contig="):
            if re.search("ID=chr([0-9]+|[X]),", line) is None:
                continue
            else:
                sys.stdout.write(line)
        else:
            sys.stdout.write(line)
    else:
        row = line.split("\t")
        chrom = row[0]
        if re.match("chr([0-9]+|[X])", chrom) is None:
            continue
        svtype = row[2].split(".")[1]
        if svtype != "DEL" and svtype != "INS":
            continue
        sys.stdout.write(line)