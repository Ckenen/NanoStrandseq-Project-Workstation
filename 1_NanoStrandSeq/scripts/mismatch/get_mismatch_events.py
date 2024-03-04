#!/usr/bin/env python
import sys
from collections import defaultdict
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import SegmentTools
from pyBioInfo.Utils import BundleBuilder


def load_segments(path, chrom):
    with pysam.AlignmentFile(path) as f:
        for segment in f.fetch(contig=chrom):
            chrom = segment.reference_name
            start = segment.reference_start
            end = segment.reference_end
            name = segment.query_name
            strand = "-" if segment.is_reverse else "+"
            obj = GRange(chrom=chrom, start=start, end=end, name=name, strand=strand)
            obj.segment = segment
            # obj.xg = segment.get_tag("XG")
            yield obj


def get_events(align):
    events = []
    for item in align.parsed_md:
        if item[0] == "X":
            mapped_start, mapped_end = item[1]
            assert mapped_end - mapped_start == 1
            ref_base = item[2]
            alt_base = None
            quality = None
            position = None
            found = False
            for d in align.parsed_cigar:
                if d[0] == "M" and d[1][0] <= mapped_start and mapped_end <= d[1][1]:
                    offset = mapped_start - d[1][0]
                    read_start, read_end = d[2]
                    chrom_start, chrom_end = d[3]
                    alt_base = align.sequence[read_start + offset]
                    quality = align.qualities[read_start + offset]
                    position = chrom_start + offset
                    events.append([position, ref_base, alt_base, quality])
                    found = True
                    break
            assert found
        elif item[0] == "D":
            mapped_start, mapped_end = item[1]
            for d in align.parsed_cigar:
                if d[0] == "D" and d[1][0] <= mapped_start and mapped_end <= d[1][1]:
                    assert d[1][0] == mapped_start
                    assert d[1][1] == mapped_end
                    for i in range(mapped_end - mapped_start):
                        ref_base = item[2][i]
                        position = d[3][0] + i
                        events.append([position, ref_base, "-", 0])
                        # print([position, ref_base, "-", 0])
                    break
            # print(item[2])
    return events


def main():
    input_bam, chrom, out_tsv = sys.argv[1:]
    loader = load_segments(input_bam, chrom)
    with open(out_tsv, "w+") as fw:
        for bundle in BundleBuilder(loader, keep=True):
            chrom = bundle.chrom
            for strand in ["+", "-"]:
                data = defaultdict(list)
                for align in bundle.data:
                    if align.strand != strand:
                        continue
                    align.sequence = align.segment.query_sequence
                    align.qualities = align.segment.query_qualities
                    align.parsed_md = SegmentTools.parse_md_tag(align.segment)
                    align.parsed_cigar = SegmentTools.parse_cigar(align.segment)
                    events = get_events(align)
                    # [position, ref_base, alt_base, quality]
                    for e in events:
                        data[e[0]].append((e[1], e[2], e[3]))
                for position, items in data.items():
                    s = ";".join([",".join(map(str, item)) for item in items])
                    line = "\t".join(map(str, [chrom, position, strand, s]))
                    fw.write(line + "\n")


if __name__ == "__main__":
    main()
