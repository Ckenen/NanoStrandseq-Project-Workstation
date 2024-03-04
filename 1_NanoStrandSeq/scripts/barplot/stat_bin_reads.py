#!/usr/bin/env python
import re
import optparse
import multiprocessing
import pandas as pd
import pysam


def stat_bin_reads(bamfile, chrom, width, rmdup, rmlc):
    with pysam.AlignmentFile(bamfile) as f:
        length = f.get_reference_length(chrom)
        nbin = int(length / width)
        if length % width > 0:
            nbin += 1
        
        rows = []
        offset = 5
        for i in range(nbin):
            start = i * width
            end = min(start + width, length)
            row = [chrom, i, start, end, end - start, 0, 0, 0, 0, 0, 0]
            rows.append(row)
            
        for s in f.fetch(chrom):
            isdup = s.is_duplicate
            if rmdup and isdup:
                continue
            
            islc = False
            if s.has_tag("XH"):
                islc = s.get_tag("XH") == "N"
            if rmlc and islc:
                continue
                
            start = s.reference_start
            end = s.reference_end
            strand = "-" if s.is_reverse else "+"
            if s.is_paired and s.is_read2: # Paired-end
                strand = "-" if strand == "+" else "+"
                
            parental = "U"
            if s.has_tag("XP"):
                parental = s.get_tag("XP")
                
            idx = int(start / width)
            if strand == "+":
                rows[idx][offset] += 1
                if parental == "P":
                    rows[idx][offset + 1] += 1
                elif parental == "M":
                    rows[idx][offset + 2] += 1
            else:
                rows[idx][offset + 3] += 1
                if parental == "P":
                    rows[idx][offset + 4] += 1
                elif parental == "M":
                    rows[idx][offset + 5] += 1
                
        columns = ["Chrom", "Bin", "Start", "End", "Width", 
               "Crick", "Crick.P", "Crick.M", 
               "Watson", "Watson.P", "Watson.M"]
        df = pd.DataFrame(rows, columns=columns)
        
        return df


def main():
    parser = optparse.OptionParser(usage="%prog [options] input.bam output.tsv")
    parser.add_option("-d", "--remove-duplicates", dest="rmdup", action="store_true", default=False, 
                      help="Remove duplicate reads.")
    parser.add_option("-c", "--remove-low-confidence", dest="rmlc", action="store_true", default=False, 
                      help="Remove low confidence reads.")
    parser.add_option("-w", "--bin-width", dest="width", type="int", default=1000000)
    parser.add_option("-t", "--threads", dest="threads", type="int", default=1)
    options, args = parser.parse_args()
    
    bamfile, tsvfile = args
    width = options.width
    threads = options.threads
    rmdup = options.rmdup
    rmlc = options.rmlc
    
    results = []
    pool = multiprocessing.Pool(threads)
    with pysam.AlignmentFile(bamfile) as f:
        for chrom in f.references:
            args = (bamfile, chrom, width, rmdup, rmlc)
            r = pool.apply_async(stat_bin_reads, args)
            results.append(r)
    pool.close()
    pool.join()
    
    dats = [r.get() for r in results]
    dat = pd.concat(dats, axis=0)
    dat.to_csv(tsvfile, sep="\t", index=False)
           
            
if __name__ == "__main__":
    main()