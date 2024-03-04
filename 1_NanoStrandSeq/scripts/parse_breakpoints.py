#!/usr/bin/env python
import sys
import re
from collections import defaultdict
import multiprocessing as mp
import pysam


def parse_cigar_string(s):
    array = []
    while len(s) > 0:
        ret = re.search("^[0-9]+[MIDSH]", s)
        i1, i2, = ret.span()
        s1 = s[:i2]
        k, v = s1[-1], int(s1[:-1])
        s = s[i2:]
        array.append([k ,v])
    return array


def get_read_length(cigars):
    v = 0
    for cigar in cigars:
        if cigar[0] in ["M", "S", "H", "I"]:
            v += cigar[1]
    return v


def get_clips(cigars):
    clip1 = 0
    clip2 = 0
    if len(cigars) >= 2:
        if cigars[0][0] in ["S", "H"]:
            clip1 = cigars[0][1]
        if cigars[-1][0] in ["S", "H"]:
            clip2 = cigars[-1][1]
    return clip1, clip2


def get_mapped_index(cigars, strand, length):
    clip1, clip2 = get_clips(cigars)
    if strand == "-":
        clip1, clip2 = clip2, clip1
    return [clip1, length - clip2] 


def get_mapped_region(cigars, start):
    end = start
    for cigar in cigars:
        if cigar[0] in ["M", "D", "N"]:
            end += cigar[1]
    return start, end
    
    
def parse_supplementary_alignment(s):
    d = dict()
    chrom, start, strand, cigarstring, mapq, nm = s.split(",")
    start, mapq, nm = int(start), int(mapq), int(nm)
    cigars = parse_cigar_string(cigarstring)
    length = get_read_length(cigars)
    start, end = get_mapped_region(cigars, start)
    i1, i2 = get_mapped_index(cigars, strand, length)
    d["Chrom"] = chrom
    d["Start"] = start
    d["End"] = end
    d["Strand"] = strand
    d["Length"] = length
    d["CigarString"] = cigarstring
    d["Cigars"] = cigars
    d["ReadStart"] = i1
    d["ReadEnd"] = i2
    d["MapQ"] = mapq
    d["NM"] = nm
    return d


def worker(f_bam, chrom, start, end):
    # print(chrom, start, end, sep="\t")
    counter = defaultdict(int)
    with pysam.AlignmentFile(f_bam) as f:
        for s in f.fetch(chrom, start, end):
            if s.reference_start < start:
                continue
            hits = []

            # primary
            cigars = parse_cigar_string(s.cigarstring)
            length = get_read_length(cigars)
            strand = "+" if s.is_forward else "-"
            i1, i2 = get_mapped_index(cigars, strand, length)
            hit1 = dict()
            hit1["Chrom"] = s.reference_name
            hit1["Start"] = s.reference_start
            hit1["End"] = s.reference_end
            hit1["Strand"] = strand
            hit1["Length"] = length
            hit1["CigarString"] = s.cigarstring
            hit1["Cigars"] = cigars
            hit1["ReadStart"] = i1
            hit1["ReadEnd"] = i2
            hit1["MapQ"] = s.mapping_quality
            hit1["NM"] = s.get_tag("NM")
            hits.append(hit1)
            
            min_clip = 200
            clip1, clip2 = get_clips(hit1["Cigars"])
            if clip1 < min_clip and clip2 < min_clip:
                continue
            if clip1 >= min_clip and clip2 >= min_clip:
                continue
            
            if True:
                if s.has_tag("SA"):
                    for sa in s.get_tag("SA").split(";"):
                        if sa == "":
                            continue
                        hit2 = parse_supplementary_alignment(sa)
                        if hit2["Chrom"] == s.reference_name:
                            hits.append(hit2)
            if len(hits) != 2:
                continue
                
            hits = list(sorted(hits, key=lambda item: item["Start"]))
            
            hit1, hit2 = hits
            clip1, clip2 = get_clips(hit1["Cigars"])
            clip3, clip4 = get_clips(hit2["Cigars"])
            if clip1 >= min_clip and clip2 < min_clip and clip3 >= min_clip and clip4 < min_clip:
                counter[hit1["Start"]] += 1
                counter[hit2["Start"]] += 1
            elif clip1 < min_clip and clip2 >= min_clip and clip3 < min_clip and clip4 >= min_clip:
                counter[hit1["End"]] += 1
                counter[hit2["End"]] += 1
            elif clip1 >= min_clip and clip2 < min_clip and clip3 < min_clip and clip4 >= min_clip:
                counter[hit1["Start"]] += 1
                counter[hit2["End"]] += 1
            elif clip1 < min_clip and clip2 >= min_clip and clip3 >= min_clip and clip4 < min_clip:
                counter[hit1["End"]] += 1
                counter[hit2["Start"]] += 1
    return [chrom, counter]


def main():
    f_bam, threads, f_out = sys.argv[1:]
    threads = int(threads)
    
    results = []
    pool = mp.Pool(threads)
    with pysam.AlignmentFile(f_bam) as f:
        chroms = list(f.references)
        for chrom in chroms:
            length = f.get_reference_length(chrom)
            step = 10000000
            for start in range(0, length, step):
                end = min(start + step, length)
                args = (f_bam, chrom, start, end)
                r = pool.apply_async(worker, args)
                results.append(r)
                # break
    pool.close()
    pool.join()
    
    counter_all = defaultdict(int)
    for r in results:
        chrom, counter = r.get()
        for k, v in counter.items():
            counter_all[(chrom, k)] += v
            
    with open(f_out, "w+") as fw:
        for (chrom, start), count in sorted(counter_all.items()):
            fw.write("%s\t%d\t%d\t%d\n" % (chrom, start, start + 1, count))
    
    
if __name__ == "__main__":
    main()
    
    