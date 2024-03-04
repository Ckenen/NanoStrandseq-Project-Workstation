#!/usr/bin/env python
import sys
import re
import pysam


CHROM_PATTERN = "^chr([0-9]+|[XY])$"
MIN_MAPQ = 30

    
def filter_sa_tag(segment, pattern):
    try:
        sa_tag_value = segment.get_tag("SA")
        items = []
        for s in sa_tag_value.split(";"):
            if s == "":
                continue
            if re.match(pattern, s.split(",")[0]) is None:
                continue
            items.append(s)
        if len(items) == 0:
            sa_tag_value = None
        else:
            sa_tag_value = ";".join(items) + ";"
        segment.set_tag("SA", sa_tag_value)
    except KeyError:
        pass
    return segment


def main():
    infile, outfile = sys.argv[1:]

    with pysam.AlignmentFile(infile) as f, \
        pysam.AlignmentFile(outfile, "wb", f) as fw:
        for chrom in f.references:
            if re.match(CHROM_PATTERN, chrom.split(",")[0]) is None:
                continue
            for segment in f.fetch(chrom):
                if segment.mapping_quality < MIN_MAPQ:
                    continue
                if segment.is_secondary:
                    continue
                # if segment.is_supplementary:
                #     continue
                segment = filter_sa_tag(segment, CHROM_PATTERN)
                fw.write(segment)


if __name__ == "__main__":
    main()