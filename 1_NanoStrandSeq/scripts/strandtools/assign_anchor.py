#!/usr/bin/env python
import sys
import json
from collections import Counter, defaultdict
import numpy as np
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import BedFile
from pyBioInfo.Utils import ShiftLoader, BundleBuilder, SegmentTools


def load_reads(path, chrom, names):
    with pysam.AlignmentFile(path) as f:
        for segment in f.fetch(chrom):
            dn = segment.get_tag("DN")
            if dn not in names:                
                continue
            start = segment.reference_start
            end = segment.reference_end
            strand = "-" if segment.is_reverse else "+"
            obj = GRange(chrom=chrom, start=start, end=end, strand=strand)
            obj.events = SegmentTools.get_events(segment)
            yield obj


def main():
    bamfile, bedfile, txtfile = sys.argv[1:]

    with open(txtfile) as f:
        duplicate_set_names = json.load(f)
    
    chrom_hets = defaultdict(list)
    with BedFile(bedfile) as f:
        for obj in f:
            obj.allele1 = obj.name[0]
            obj.allele2 = obj.name[2]
            chrom_hets[obj.chrom].append(obj)
                
    for chrom, names in duplicate_set_names.items():
        hets = chrom_hets[chrom]
        if len(hets) == 0:
            continue
            
        names = set(names)
        if len(names) < 1:
            continue
        
        loader = ShiftLoader(hets)
        reads = load_reads(bamfile, chrom, names)
        for bundle in BundleBuilder(reads, keep=True):
            if len(bundle.data) < 2: # 3
                continue
            chrom = bundle.chrom
            start = bundle.start_min
            end = bundle.end_max
            reads_crick = [] # crick
            reads_watson = [] # watson
            hets_local = list(loader.fetch(chrom=chrom, start=start, end=end))
            if len(hets_local) == 0:
                continue
            for obj in bundle.data:
                if obj.strand == "+":
                    reads_crick.append(obj)
                else:
                    reads_watson.append(obj)
            for reads, strand in zip([reads_crick, reads_watson], ["+", "-"]):
                if len(reads) < 2: # 3
                    continue
                coverages = np.zeros(end - start, dtype=np.int)
                meta = defaultdict(list)
                for obj in reads:
                    for idx in range(obj.start - start, obj.end - start):
                        coverages[idx] += 1
                    for e in obj.events:
                        idx = e[0] - start
                        if e[1] == "-":
                            continue
                        for bi in range(len(e[1])):
                            meta[e[0] + bi].append(e[2])
                for anchor in hets_local:
                    cov = coverages[anchor.start - start]
                    if cov < 2: # 3
                        continue
                    ref = anchor.name[0]
                    c = Counter(meta[anchor.start])
                    alt_count = sum(c.values()) # altanative base count
                    del_count = c["-"]
                    ref_count = cov - alt_count
                    base_count = cov - del_count
                    if base_count < 1:
                        continue
                    cutoff = int(base_count * 0.75) # 0.8
                    if cutoff < base_count * 0.75:
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


if __name__ == "__main__":
    main()