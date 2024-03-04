#!/usr/bin/env python
import sys
import numpy as np
import pysam


def load_coverage(path):
    data = dict()
    wbin = 1000000
    with pysam.AlignmentFile(path) as f:
        for item in f.header.as_dict()["SQ"]:
            sn = item["SN"]
            ln = item["LN"]
            nbin = int(ln / wbin)
            if ln % wbin > 0:
                nbin += 1
            data[sn] = [np.zeros(nbin), np.zeros(nbin)]
        for segment in f:
            if segment.get_tag("XT") == "F":
                if segment.is_reverse:
                    s = 1
                else:
                    s = 0
            else:
                if segment.is_reverse:
                    s = 0
                else:
                    s = 1
            counts = data[segment.reference_name][s]
            a = int(segment.reference_start / wbin)
            b = int((segment.reference_end - 1) / wbin)
            c = b - a + 1
            w = 1 / c
            for x in range(a, b + 1):
                counts[x] += w
    return data


def main():
    infile, outfile = sys.argv[1:]

    data = load_coverage(infile)

    with open(outfile, "w+") as fw:
        for chrom in sorted(data):
            values1, values2 = data[chrom]
            line1 = "%s\t+\t%s\n" % (chrom, ",".join(map(str, values1)))
            line2 = "%s\t-\t%s\n" % (chrom, ",".join(map(str, values2)))
            fw.write(line1)
            fw.write(line2)


if __name__ == "__main__":
    main()
