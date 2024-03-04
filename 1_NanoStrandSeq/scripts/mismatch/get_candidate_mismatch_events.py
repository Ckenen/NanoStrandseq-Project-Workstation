#!/usr/bin/env python
import sys
import gzip
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import ShiftLoader


def load_mismatch_events(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            chrom, position, xg, events = line.strip("\n").split("\t")
            position = int(position)
            obj = GRange(chrom=chrom, start=position, end=position + 1)
            obj.xg = xg
            obj.events = events
            yield obj


def main():
    in_tsv1, in_tsv2, out_tsv = sys.argv[1:]

    loader = ShiftLoader(load_mismatch_events(in_tsv2))
    with gzip.open(in_tsv1, "rt") as f, gzip.open(out_tsv, "wt") as fw:
        for line in f:
            chrom, position = line.strip("\n").split("\t")
            start = int(position)
            events1 = ""
            events2 = ""
            for obj in loader.fetch(chrom=chrom, start=start, end=start + 1):
                if obj.xg == "C":
                    events1 = obj.events
                else:
                    events2 = obj.events
            fw.write("%s|%s\n" % (events1, events2))


if __name__ == "__main__":
    main()
