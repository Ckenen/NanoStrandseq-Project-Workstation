#!/usr/bin/env python
import sys
import numpy as np
from collections import defaultdict
from pyBioInfo.IO.File import BamFile
from pyBioInfo.Utils import SegmentTools, BundleBuilder


def load_alignments(path):
    with BamFile(path) as bam:
        for alignment in bam:
            if alignment.segment.get_tag("XH") == "Y":
                yield alignment


def get_coverage(bundle):
    width = bundle.end_max - bundle.start_min
    coverage = np.zeros(width, dtype=np.int)
    for alignment in bundle.data:
        for cigar in SegmentTools.parse_cigar(alignment.segment):
            if cigar[0] == "M":
                for i in range(cigar[3][0] - bundle.start_min, cigar[3][1] - bundle.start_min):
                    coverage[i] += 1
    return coverage


def get_mismatch_events(bundle):
    events = []
    for alignment in bundle.data:
        segment = alignment.segment
        parsed_cigar = SegmentTools.parse_cigar(segment)
        parsed_md = SegmentTools.parse_md_tag(segment)
        sequence = segment.query_sequence
        qualities = segment.query_qualities
        for md in parsed_md:
            if md[0] == "X":
                ref_base = md[2]
                hit = False
                for cigar in parsed_cigar:
                    if cigar[0] == "M" and cigar[1][0] <= md[1][0] and cigar[1][1] >= md[1][1]:
                        offset = md[1][0] - cigar[1][0]
                        read_offset = offset + cigar[2][0]
                        position = offset + cigar[3][0]
                        alt_base = sequence[read_offset]
                        quality = qualities[read_offset]
                        hit = True
                        break
                assert hit
                events.append((position, ref_base, alt_base, quality))
    return events


def process_bundle(bundle, fw):
    coverage = get_coverage(bundle)
    events = get_mismatch_events(bundle)
    data = defaultdict(list)
    for event in events:
        data[event[0]].append(event)
    for position in sorted(data.keys()):
        tmp = data[position]
        cov = coverage[position - bundle.start_min]
        num = 0
        num = len(tmp)
        ratio = num / cov
        assert ratio >= 0
        assert ratio <= 1
        if ratio >= 0.2:
            line = "\t".join(
                map(str, [bundle.chrom, position, position + 1, "%d/%d" % (num, cov), ratio, "+"]))
            fw.write(line + "\n")


def main():
    in_bam, out_bed = sys.argv[1:]

    loader = load_alignments(in_bam)
    with open(out_bed, "w+") as fw:
        for bundle in BundleBuilder(loader, min_capacity=0, keep=True):
            if len(bundle.data) < 10:
                continue
            # print(bundle.as_dict())
            # print(len(bundle.data))
            process_bundle(bundle, fw)


if __name__ == '__main__':
    main()
