#!/usr/bin/env python
import sys
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import BedFile
from pyBioInfo.Utils import ShiftLoader


def load_hets_from_vcf(path):
    hets = []
    with pysam.VariantFile(path) as f:
        sample = f.header.samples[0]
        for record in f:
            gt = record.samples[sample]["GT"]
            allele1 = record.alleles[gt[0]]
            allele2 = record.alleles[gt[1]]
            if len(allele1) > 1 or len(allele2) > 1 or allele1 == allele2:
                continue
            start, end = record.pos - 1, record.pos
            start = end - 1
            het = GRange(chrom=record.chrom, start=start, end=end, name="HET")
            het.allele1 = allele1
            het.allele2 = allele2
            hets.append(het)
    hets.sort()
    return hets


def filter_hets(hets, regions):
    loader = ShiftLoader(regions)
    hets1 = []
    for het in hets:
        n = len(list(loader.fetch(obj=het)))
        if n > 0:
            hets1.append(het)
    return hets1


def get_snv_set(hets):
    array = []
    for het in hets:
        a1, a2 = het.allele1, het.allele2
        if a1 > a2:
            a1, a2 = a2, a1
        array.append((het.chrom, het.start, a1, a2))
    return set(array)


def main():
    # ref.vf query.vcf regions.bed output.txt
    vcffile1, vcffile2, bedfile, outfile = sys.argv[1:]

    regions = []
    with BedFile(bedfile) as f:
        regions = [x for x in f]
    hets1 = load_hets_from_vcf(vcffile1) # reference
    hets2 = load_hets_from_vcf(vcffile2) # query
    hets1 = filter_hets(hets1, regions)
    hets2 = filter_hets(hets2, regions)
    hets1 = get_snv_set(hets1)
    hets2 = get_snv_set(hets2)

    hets3 = hets1 & hets2
    precision = len(hets3) / len(hets2)
    recall = len(hets3) / len(hets1)
    f1 = 2 / (1 / precision + 1 / recall)
    with open(outfile, "w+") as fw:
        fw.write("Hets1\tHets2\tHetsCommon\tPrecision\tRecall\tF1\n")
        fw.write("%d\t%d\t%d\t%f\t%f\t%f\n" % (len(hets1), len(hets2), len(hets3), precision, recall, f1))

    with open(outfile + ".references", "w+") as fw:
        for snv in hets1:
            fw.write("\t".join(map(str, snv)) + "\n")
    with open(outfile + ".queries", "w+") as fw:
        for snv in hets2:
            fw.write("\t".join(map(str, snv)) + "\n")
    with open(outfile + ".commons", "w+") as fw:
        for snv in hets3:
            fw.write("\t".join(map(str, snv)) + "\n")


if __name__ == "__main__":
    main()