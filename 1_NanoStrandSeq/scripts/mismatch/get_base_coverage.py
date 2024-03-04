#!/usr/bin/env python
import sys
import numpy as np
from pyBioInfo.IO.File import BamFile
from pyBioInfo.Utils import SegmentTools, BundleBuilder


def load_alignments(path):
    with BamFile(path) as f:
        for align in f:
            xg = align.segment.get_tag("XG")
            align.xg = xg
            yield align


def get_coverage(bundle):
    bundle_start = bundle.start_min
    length = bundle.end_max - bundle.start_min
    covs1 = np.zeros(length, dtype=np.int)  # Crick
    covs2 = np.zeros(length, dtype=np.int)  # Watson
    covs = None
    for align in bundle.data:
        if align.xg == "Crick":
            covs = covs1
        else:
            covs = covs2
        parsed_cigar = SegmentTools.parse_cigar(align.segment)

        blocks = []
        for item in parsed_cigar:
            if item[0] == "M":
                chrom_start, chrom_end = item[3]
                blocks.append([chrom_start, chrom_end])

        for i1, i2 in blocks:
            for i in range(i1 - bundle_start, i2 - bundle_start):
                covs[i] += 1
    return covs1, covs2


def collapse(values):
    array = []
    start = -1
    last_v = 0
    for i, v in enumerate(values):
        if v == 0:
            if last_v == 0:
                continue
            else:
                array.append([start, i, last_v])
        else:
            if last_v == 0:
                start = i
            elif v == last_v:
                continue
            else:
                array.append([start, i, last_v])
                start = i
        last_v = v
    if last_v != 0:
        array.append([start, i + 1, last_v])
    return array


def main():
    input_bam, out_c_bg, out_w_bg = sys.argv[1:]
    loader = load_alignments(input_bam)
    with open(out_c_bg, "w+") as fw1, open(out_w_bg, "w+") as fw2:
        for bundle in BundleBuilder(loader, min_capacity=0, keep=True):
            covs1, covs2 = get_coverage(bundle)

            for i1, i2, v in collapse(covs1):
                p1, p2 = bundle.start_min + i1, bundle.start_min + i2
                line = "\t".join(map(str, [bundle.chrom, p1, p2, v]))
                fw1.write(line + "\n")
            for i1, i2, v in collapse(covs2):
                p1, p2 = bundle.start_min + i1, bundle.start_min + i2
                line = "\t".join(map(str, [bundle.chrom, p1, p2, v]))
                fw2.write(line + "\n")


if __name__ == "__main__":
    main()
