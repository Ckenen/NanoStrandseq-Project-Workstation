#!/usr/bin/env python
import sys, re
import numpy as np
import pysam

PATTERN = "^chr[0-9]+$"

def infer_library_layout(path):
    is_paired = False
    with pysam.AlignmentFile(path) as f:
        for segment in f:
            if segment.is_unmapped:
                continue
            is_paired = segment.is_paired
            break
    return is_paired
    

def load_fragments(path):
    is_paired = infer_library_layout(path)
    with pysam.AlignmentFile(path) as f:
        for chrom in f.references:
            if re.match(PATTERN, chrom) is None:
                continue
            segments = []
            for segment in f.fetch(chrom):
                if segment.is_unmapped:
                    continue
                if is_paired and (not segment.is_proper_pair):
                    continue
                if segment.is_secondary:
                    continue
                if segment.is_supplementary:
                    continue
                segments.append(segment)                     
            if is_paired:
                segments = list(sorted(segments, key=lambda item: item.query_name))
                i = 0
                while i < len(segments) - 1:
                    segment1 = segments[i]
                    segment2 = segments[i + 1]
                    if segment1.query_name == segment2.query_name:
                        assert segment1.next_reference_start == segment2.reference_start
                        assert segment2.next_reference_start == segment1.reference_start
                        start = min(segment1.reference_start, segment2.reference_start)
                        end = max(segment1.reference_end, segment2.reference_end)
                        yield [chrom, start, end, segment1, segment2]                        
                        i += 2
                    else:
                        i += 1
            else:
                for segment in segments:
                    start = segment.reference_start
                    end = segment.reference_end
                    yield [chrom, start, end, segment]
    
    
def main():
    infile, outfile = sys.argv[1:]
    
    assert infile.endswith(".bam")
    cell = infile.split("/")[-1][:-4]
    
    lengths = []
    for items in load_fragments(infile):
        start, end = items[1:3]
        lengths.append(end - start)
        
    with open(outfile, "w+") as fw:
        for length in lengths:
            fw.write("%d\n" % length)
                
    mean = 0
    median = 0
    std = 0
    if len(lengths) > 0:
        mean = np.mean(lengths)
        median = np.median(lengths)
        std = np.std(lengths)
    print("Cell\tCount\tMean\tMedian\tStd")
    print(cell, len(lengths), mean, median, std, sep="\t")

    
if __name__ == '__main__':
    main()
