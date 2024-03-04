#!/usr/bin/env python
import sys, os, pysam
import pandas as pd
from pyBioInfo.IO.File import BedFile
from pyBioInfo.Utils import ShiftLoader


def load_snvs(path):
    snvs = dict()
    with pysam.VariantFile(path) as f:
        sample = list(f.header.samples)[0]
        for record in f:
            gt = record.samples[sample]["GT"]
            ps = ""
            try:
                ps = record.samples[sample]["PS"]
            except KeyError:
                pass
            a1, a2 = record.alleles[gt[0]], record.alleles[gt[1]]
            if len(a1) > 1 or len(a2) > 1:
                continue
            snvs[(record.chrom, record.start)] = [a1, a2, ps]
    return snvs

    
def filter_snvs(snvs, regions):
    new_snvs = dict()
    loader = ShiftLoader(regions)
    for (chrom, start), v in sorted(snvs.items()):
        if len(list(loader.fetch(chrom=chrom, start=start, end=start + 1))) > 0:
            new_snvs[(chrom, start)] = v
    return new_snvs


def main():
    vcf1, vcf2, bed, outdir = sys.argv[1:]
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    # Reference SNVs
    snvs1 = load_snvs(vcf1)
    # Query SNVs
    snvs2 = load_snvs(vcf2)
    # Benchmark regions
    if bed == "null":
        regions = None
    else:
        with BedFile(bed) as f:
            regions = [x for x in f]
        
    print("Loaded %d SNVs in reference VCF." % len(snvs1))
    print("Loaded %d SNVs in query VCF." % len(snvs2))
    
    if regions is not None:
        snvs1 = filter_snvs(snvs1, regions)
        snvs2 = filter_snvs(snvs2, regions)
    print("After filter, there is %d SNVs in reference." % len(snvs1))
    print("After filter, there is %d SNVs in query." % len(snvs2))
    
    inconsistents = []
    
    # Benchmark of calling
    sites = list(sorted(snvs1.keys() & snvs2.keys()))
    n1 = len(snvs1)
    n2 = len(snvs2)
    n3 = len(sites)
    r = n3 / n1
    p = n3 / n2
    print("Benchmark of calling:")
    print("Reference: %d, query: %d, overlap: %d, recall: %f, preicision: %f" % (n1, n2, n3, r, p))
    with open(outdir + "/benchmark_of_calling.tsv", "w+") as fw:
        fw.write("Reference\tQuery\tOverlap\tRecall\tPrecision\n")
        fw.write("\t".join(map(str, [n1, n2, n3, r, p])) + "\n")
    for chrom, start in (snvs1.keys() - snvs2.keys()):
        inconsistents.append([chrom, start, "Reference.only"])
    for chrom, start in (snvs2.keys() - snvs1.keys()):
        inconsistents.append([chrom, start, "Query.only"])
    
    # Benchmark of genotyping
    rows = [
        ["HOM", "HOM", 0],
        ["HOM", "HOM.2", 0],
        ["HOM", "HET", 0],
        ["HET", "HOM", 0],
        ["HET", "HET", 0],
        ["HET", "HET.2", 0]]
    for chrom, start in sites:
        snv1 = snvs1[(chrom, start)]
        snv2 = snvs2[(chrom, start)]
        n1 += 1
        a1, a2 = snv1[0], snv1[1]
        b1, b2 = snv2[0], snv2[1]
        if a1 == a2:
            if b1 == b2:
                if a1 == b1:
                    rows[0][2] += 1 # HOM HOM
                else:
                    rows[1][2] += 1
                    inconsistents.append([chrom, start, "Diff.genotype"])
            else:
                rows[2][2] += 1
                inconsistents.append([chrom, start, "Diff.genotype"])
        else: # a1 != a2
            if b1 == b2:
                rows[3][2] += 1
                inconsistents.append([chrom, start, "Diff.genotype"])
            else: # b1 != b2
                if (a1 == b1 and a2 == b2) or (a1 == b2 and a2 == b1):
                    rows[4][2] += 1 # HET HET
                else:
                    rows[5][2] += 1
                    inconsistents.append([chrom, start, "Diff.genotype"])
    vs = [row[2] for row in rows]
    n1 = sum(vs)
    n2 = vs[0] + vs[4]
    p = n2 / n1
    print("Benchmark of genotyping:")
    print("Total: %d, identical: %d, precision: %f" % (n1, n2, p))
    with open(outdir + "/benchmark_of_genotyping.tsv", "w+") as fw:
        fw.write("Total\tIdentical\tPrecision\n")
        fw.write("\t".join(map(str, [n1, n2, p])) + "\n")
    df = pd.DataFrame(rows)
    df.columns = ["Genotype1", "Genotype2", "Number"]
    df.to_csv(outdir + "/benchmark_of_genotyping.detail.tsv", sep="\t", index=False)    
    
    # Benchmark of phasing
    n1 = 0
    n2 = 0
    tmp1 = []
    tmp2 = []
    for chrom, start in sites:
        snv1 = snvs1[(chrom, start)]
        snv2 = snvs2[(chrom, start)]
        if snv1[0] == snv1[1]:
            continue
        if snv2[0] == snv2[1]:
            continue
        if snv1[2] != "PATMAT":
            continue
        if snv2[2] != "PATMAT":
            continue
        if snv1[0] == snv2[0] and snv1[1] == snv2[1]:
            n1 += 1
            tmp1.append([chrom, start, "Diff.phase"])
        else:
            n2 += 1
            tmp2.append([chrom, start, "Diff.phase"])
    tmp = []
    if n1 > n2:
        tmp = tmp2
    elif n1 < n2:
        tmp = tmp1
    inconsistents.extend(tmp)
    p = 0
    if n1 + n2 > 0:
        p = max(n1, n2) / (n1 + n2)
    print("Benchmark of phasing:")
    print("Total: %d, identical: %d, precision: %f" % (n1 + n2, max(n1, n2), p))
    with open(outdir + "/benchmark_of_phasing.tsv", "w+") as fw:
        fw.write("Total\tHP1_HP1\tHP1_HP2\tIdentical\tPrecision\n")
        fw.write("\t".join(map(str, [n1 + n2, n1, n2, max(n1, n2), p])) + "\n")
    
    # Inconsistents
    path = outdir + "/inconsistents.bed"
    with open(path, "w+") as fw:
        for chrom, start, name in sorted(inconsistents):
            color = "255,0,0"
            if name == "Reference.only":
                color = "0,255,0"
            if name == "Query.only":
                color = "0,0,255"
            if name == "Diff.genotype":
                color = "0,0,0"
            line = "\t".join(map(str, [chrom, start, start + 1, name, ".", "+", start, start + 1, color]))
            fw.write(line + "\n")
    assert os.system("bgzip %s" % path) == 0
    assert os.system("tabix -p bed %s.gz" % path) == 0
    
    
if __name__ == "__main__":
    main()
    