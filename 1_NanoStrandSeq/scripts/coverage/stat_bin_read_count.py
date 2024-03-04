#!/usr/bin/env python
import optparse
import numpy as np
import pysam

usage = """
    %prog [options] input.bam output.txt
"""


def main():
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-w", "--width", dest="width",
                      type="int", default=200000, help="")
    parser.add_option("-c", "--high-conf", dest="high_conf",
                      action="store_true", default=False, help="Only consider high-confidence alignments. (XH:Z:Y)")
    parser.add_option("-d", "--rmdup", dest="rmdup",
                      action="store_true", default=False, help="Only consider not pcr duplicated alignments. (DT:Z:N)")
    parser.add_option("-p", "--phased", dest="phased", 
                      action="store_true", default=False, help="Output paternal and maternal read count. (XP:Z:P, XP:Z:M)")
    (option, argv) = parser.parse_args()
    if len(argv) != 2:
        parser.print_help()
        exit(1)

    infile, outfile = argv
    bin_width = option.width
    assert bin_width > 0
    high_conf = option.high_conf
    rmdup = option.rmdup
    phased = option.phased

    with pysam.AlignmentFile(infile) as f:
        data = dict()
        for chrom, length in zip(f.references, f.lengths):
            c = int(length / bin_width)
            if length % bin_width > 0:
                c += 1
            # crick watson crick.p crick.m watson.p watson.m
            data[chrom] = [np.zeros(c), np.zeros(c), 
                           np.zeros(c), np.zeros(c), 
                           np.zeros(c), np.zeros(c)]
        for segment in f:
            chrom = segment.reference_name
            i1 = int(segment.reference_start / bin_width)
            i2 = int((segment.reference_end - 1) / bin_width)
            v = 1 / (i2 - i1 + 1)

            source = "U"
            if phased:
                try:
                    source = segment.get_tag("XP")
                except KeyError:
                    source = "U"
            if high_conf:
                try:
                    xh = segment.get_tag("XH")
                    if xh != "Y":
                        continue
                except KeyError:
                    pass
            if rmdup:
                if segment.get_tag("DT") != "N":
                    continue
          
            for i in range(i1, i2 + 1):
                if not segment.is_reverse:
                    data[chrom][0][i] += v
                    if source == "P":
                        data[chrom][2][i] += v
                    elif source == "M":
                        data[chrom][3][i] += v
                else:
                    data[chrom][1][i] += v
                    if source == "P":
                        data[chrom][4][i] += v
                    elif source == "M":
                        data[chrom][5][i] += v

    with open(outfile, "w+") as fw:
        for chrom in sorted(data.keys()):
            values = data[chrom]
            if not phased:
                values = values[:2]
            s = [",".join(map(str, vs)) for vs in values]
            s = "\t".join(s)
            line = "\t".join([chrom, str(bin_width), s])
            fw.write(line + "\n")


if __name__ == '__main__':
    main()
