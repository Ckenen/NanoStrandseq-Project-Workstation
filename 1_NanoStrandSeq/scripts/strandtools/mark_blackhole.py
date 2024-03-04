#!/usr/bin/env python
import sys
import pysam
from pyBioInfo.Range import GRange
# from pyBioInfo.IO.File import BamFile, BedFile
from pyBioInfo.Utils import ShiftLoader


def load_blackholes(bf, chrom):
    try:
        for line in bf.fetch(chrom):
            chrom, start, end = line.strip("\n").split("\t")[:3]
            start = int(start)
            end = int(end)
            yield GRange(chrom=chrom, start=start, end=end)
    except ValueError:
        pass
    
    
def main():
    # input.bam blackholes.bed.gz output.bam
    infile1, infile2, outfile = sys.argv[1:]
    
    n1, n2 = 0, 0
    bf = pysam.TabixFile(infile2)
    with pysam.AlignmentFile(infile1) as f, pysam.AlignmentFile(outfile, "wb", f) as fw:
        for chrom in f.references:
            loader = ShiftLoader(load_blackholes(bf, chrom))
            for segment in f.fetch(chrom):
                start = segment.reference_start
                end = segment.reference_end
                hits = list(loader.fetch(chrom=chrom, start=start, end=end))
                n1 += 1
                if len(hits) > 0:
                    bh = "Y"
                    n2 += 1
                else:
                    bh = "N"
                segment.set_tag("BH", bh)
                fw.write(segment)
    bf.close()
    
    print("%d reads" % n1)
    print("%d (%.2f%%) hit blackhole" % (n2, n2 * 100 / n1))
    
    
if __name__ == '__main__':
    main()
    