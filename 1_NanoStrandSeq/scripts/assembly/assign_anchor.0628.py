#!/usr/bin/env python
import sys
from collections import Counter, defaultdict
import numpy as np
import pysam
from pyBioInfo.IO.File import BedFile, Alignment
from pyBioInfo.Utils import ShiftLoader, BundleBuilder


def load_alignments(path, chrom):
    with pysam.AlignmentFile(path) as f:
        for segment in f.fetch(chrom):
            obj = Alignment(segment)
            events = []
            for item in segment.get_tag("ME").split(";"):
                if item == "":
                    continue
                e = item.split(",")
                e[0] = int(e[0])
                events.append(e)
            obj.events = events
            yield obj


def main():
    infile1, infile2, chrom = sys.argv[1:]
        
    with BedFile(infile2) as f:
        loader = ShiftLoader(f)
        for bundle in BundleBuilder(load_alignments(infile1, chrom), keep=True):
            if len(bundle.data) < 3:
                continue
            start = bundle.start_min
            end = bundle.end_max
            array1 = [] # crick
            array2 = [] # watson
            array3 = list(loader.fetch(chrom=chrom, start=start, end=end))
            if len(array3) == 0:
                continue
            for obj in bundle.data:
                if obj.strand == "+":
                    array1.append(obj)
                else:
                    array2.append(obj)
                    
            for array, strand in zip([array1, array2], ["+", "-"]):
                if len(array) < 3:
                    continue
                coverages = np.zeros(end - start, dtype=np.int)
                meta = defaultdict(list)
                for obj in array:
                    for idx in range(obj.start - start, obj.end - start):
                        coverages[idx] += 1
                    for e in obj.events:
                        idx = e[0] - start
                        if e[1] == "-":
                            continue
                        for bi in range(len(e[1])):
                            meta[e[0] + bi].append(e[2])
                for anchor in array3:
                    cov = coverages[anchor.start - start]
                    if cov < 3:
                        continue
                    ref = anchor.name[0]
                    c = Counter(meta[anchor.start])
                    alt_count = sum(c.values()) # altanative base count
                    del_count = c["-"]
                    ref_count = cov - alt_count
                    base_count = cov - del_count
                    cutoff = int(base_count * 0.8)
                    if cutoff < base_count * 0.8:
                        cutoff += 1
                    base = None
                    count = None
                    if ref_count >= cutoff:
                        base = ref
                        count = ref_count
                    else:
                        for k, v in c.items():
                            if k == "-":
                                continue
                            if v >= cutoff:
                                base = k
                                count = v
                                break
                    if base is None:
                        continue
                    print(chrom, anchor.start, anchor.end, "%s-%s" % (ref, base), count, strand, sep="\t")
                        

if __name__ == '__main__':
    main()
    