#!/usr/bin/env python
import sys
import pysam


def main():
    f_bam, f_out = sys.argv[1:]
    
    width = 1000000
    
    counter = dict()
    with pysam.AlignmentFile(f_bam) as f:
        for seqname in f.references:
            length = f.get_reference_length(seqname)
            nbin = int(length / width)
            if length % width > 0:
                nbin += 1
            for i in range(nbin):
                start = i * width
                end = min(start + width, length)
                counter[(seqname, i)] = [start, end, 0, 0]
                
        for s in f:
            if s.is_secondary:
                continue
            if s.is_supplementary:
                continue
            if s.mapping_quality < 30:
                continue
            if s.is_duplicate:
                continue
            strand = "+" if s.is_forward else "-"
            if s.is_paired and s.is_read2:
                strand = "+" if strand == "-" else "-"
            p = int((s.reference_start + s.reference_end) / 2)
            i = int(p / width)
            j = 2 if strand == "+" else 3
            counter[(s.reference_name, i)][j] += 1
        
    with open(f_out, "w+") as fw:
        fw.write("Name\tBin\tStart\tEnd\tCrick\tWatson\n")
        for k, v in sorted(counter.items()):
            seqname, ibin = k
            start, end, crick, watson = v
            line = "\t".join(map(str, [seqname, ibin, start, end, crick, watson]))
            fw.write(line + "\n")
    
    
if __name__ == "__main__":
    main()