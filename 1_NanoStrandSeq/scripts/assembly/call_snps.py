#!/usr/bin/env python
import sys
import pysam
import gzip
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import BedFile
from pyBioInfo.Utils import ShiftLoader
from sstools.utils import BaseMatrix2


def main():
    # in.matrix2.gz ref.vcf.gz inversion.bed.gz genome.sizes chr1 out.vcf
    infile1, infile2, infile3, infile4, chrom, outfile = sys.argv[1:]
        
    snvs1 = dict() # nanocaller
    with pysam.VariantFile(infile2) as f:
        sample = list(f.header.samples)[0]
        for record in f.fetch(chrom):
            gt = record.samples[sample]["GT"]
            a1 = record.alleles[gt[0]]
            a2 = record.alleles[gt[1]]
            if len(a1) > 1 or len(a2) > 1:
                continue
            snvs1[record.start] = [a1, a2]
    # len(snvs1)    
          
    regions = []
    with pysam.TabixFile(infile3) as f:
        for line in f.fetch(chrom):
            row = line.strip("\n").split("\t")
            chrom, start, end, name = row[:4]
            start = int(start)
            end = int(end)
            crick, watson, ratio = name.split(";")
            crick = int(crick)
            watson = int(watson)
            ratio = float(ratio)
            ps = "ROI"
            if crick >= 20 and watson >= 20 and ratio >= 0.4 and ratio <= 0.6:
                ps = "HET"
            elif watson > 40 and ratio >= 0.9:
                ps = "HOM"
            obj = GRange(chrom=chrom, start=start, end=end, name=ps)
            regions.append(obj)
    # print(len(regions))

    loader = ShiftLoader(regions)
    data = []
    with gzip.open(infile1, "rt") as f:
        for line in f:
            row = BaseMatrix2.parse_line(line)
            start, ref = row[0], row[1]
            a1, a2 = row[5], row[6]
            detail1, detail2 = row[7], row[8]
            case1 = BaseMatrix2.get_case(detail1)
            case2 = BaseMatrix2.get_case(detail2)
            conf1 = BaseMatrix2.get_confidence(a1, detail1, case1)
            conf2 = BaseMatrix2.get_confidence(a2, detail2, case2)
            row.extend([case1, conf1, case2, conf2])
            ps = "PATMAT"
            for r in loader.fetch(chrom=chrom, start=row[0], end=row[0] + 1):
                assert ps == "PATMAT"
                ps = r.name
            if case1 <= 3 and conf1:
                if case2 <= 3 and conf2:
                    data.append([chrom, start, ref, a1, a2, ps])
                else:
                    if start in snvs1:
                        b1, b2 = snvs1[start]
                        if a1 == b1:
                            a2 = b2
                        if a1 == b2:
                            a2 = b1
                        data.append([chrom, start, ref, a1, a2, ps])
            else:
                if case2 <= 3 and conf2:
                    if start in snvs1:
                        b1, b2 = snvs1[start]
                        if a2 == b1:
                            a1 = b2
                        if a2 == b2:
                            a1 = b1
                        data.append([chrom, start, ref, a1, a2, ps])
                
    with open(outfile, "w+") as fw:
        fw.write('##fileformat=VCFv4.2\n')
        with open(infile4) as f:
            for line in f:
                v1, v2 = line.strip("\n").split("\t")
                fw.write('##contig=<ID=%s,length=%s>\n' % (v1, v2))
        fw.write('##FORMAT=<ID=PS,Number=1,Type=String,Description="Phase set for GT">\n')
        fw.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Consensus Genotype across all datasets with called genotype">\n')
        fw.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNSS\n')
        for chrom, start, ref, a1, a2, ps in sorted(data):
            if ps == "HOM":
                ps = "PATMAT"
                a1, a2 = a2, a1
            alleles = [ref]
            if a1 not in alleles:
                alleles.append(a1)
            if a2 not in alleles:
                alleles.append(a2)
            gt = "%s|%s" % (alleles.index(a1), alleles.index(a2))
            line = "\t".join(map(str, [chrom, start + 1, ".", ref, ",".join(alleles[1:]), ".", "PASS", ".", "GT:PS", "%s:%s" % (gt, ps)]))
            fw.write(line + "\n")
            
            
if __name__ == "__main__":
    main()