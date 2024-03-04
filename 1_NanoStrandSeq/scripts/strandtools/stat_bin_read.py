#!/usr/bin/env python
import sys
import re
import multiprocessing
import pysam


def stat_bin_read_for_chrom(bamfile, chrom, width):
    with pysam.AlignmentFile(bamfile) as bam:
        length = bam.get_reference_length(chrom)
        nbin = int(length / width)
        if length % width > 0:
            nbin += 1
            
        data = []
        for i in range(nbin):
            start = i * width
            end = min(start + width, length)
            row = [chrom, i, start, end, end - start]
            row.extend([0] * 6 * 4)
            data.append(row)
        
        reads = []
        for segment in bam.fetch(chrom):
            start = segment.reference_start
            idx = int(start / width)
            is_dup = segment.is_duplicate
            strand = "-" if segment.is_reverse else "+"
            if segment.is_paired and segment.is_read2: # Paired-end
                strand = "-" if strand == "+" else "+"
            is_hc = True
            if segment.has_tag("XH"):
                is_hc = segment.get_tag("XH") == "Y"
            parental = "U"
            if segment.has_tag("XP"):
                parental = segment.get_tag("XP")
            reads.append([idx, strand, is_dup, is_hc, parental])
        
        conditions = [
            [False, False, 5 + 6 * 0], # All
            [True, False, 5 + 6 * 1], # RmDup
            [False, True, 5 + 6 * 2], # IsHC
            [True, True, 5 + 6 * 3]] # RmDup.IsHC
        
        for remove_dup, require_hc, offset in conditions:
            for idx, strand, is_dup, is_hc, parental in reads:
                if remove_dup and is_dup:
                    continue
                if require_hc and not is_hc:
                    continue
                if strand == "+":
                    data[idx][offset] += 1
                    if parental == "P":
                        data[idx][offset + 1] += 1
                    elif parental == "M":
                        data[idx][offset + 2] += 1
                else:
                    data[idx][offset + 3] += 1
                    if parental == "P":
                        data[idx][offset + 4] += 1
                    elif parental == "M":
                        data[idx][offset + 5] += 1
        return data
        
        
def main():
    bamfile, width, threads, outfile = sys.argv[1:]
    width = int(width)
    threads = int(threads)

    pool = None
    if threads > 1:
        pool = multiprocessing.Pool(threads)
    results = []
    with pysam.AlignmentFile(bamfile) as bam:
        for chrom in bam.references:
            if re.match("^chr([0-9]+|[XY])$", chrom) is None:
                continue
            args = (bamfile, chrom, 1000000)
            if pool:
                r = pool.apply_async(stat_bin_read_for_chrom, args)
            else:
                r = stat_bin_read_for_chrom(*args)
            results.append(r)
    if pool:
        pool.close()
        pool.join()
        results = [r.get() for r in results]
    
    columns = ["Chrom", "Bin", "Start", "End", "Width"]
    columns1 = ["Crick", "Crick.P", "Crick.M", "Watson", "Watson.P", "Watson.M"]
    columns.extend(columns1)
    columns.extend(["%s[RmDup]" % c for c in columns1])
    columns.extend(["%s[IsHC]" % c for c in columns1])
    columns.extend(["%s[RmDup.IsHC]" % c for c in columns1])
    with open(outfile, "w+") as fw:
        fw.write("\t".join(columns) + "\n")
        for rows in results:
            for row in rows:
                fw.write("\t".join(map(str, row)) + "\n")
                
                
if __name__ == "__main__":
    main()