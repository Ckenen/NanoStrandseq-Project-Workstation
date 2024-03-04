#!/usr/bin/env python
import sys
import json
import pysam
        
        
def main():
    # config.json chr1 chr1.hp1.bam chr2.hp2.bam
    txtfile, chrom, outfile1, outfile2 = sys.argv[1:]
    
    data = json.load(open(txtfile))
    
    fw1, fw2 = None, None
    cells = data["Cells"]

    for cell in cells:
        run = cell.split(".")[0]
        bamfile = "../1_NanoStrandseq/results/mapping/mark_parental/%s/%s.bam" % (run, cell)
        jsonfile = "results/%s/round2/splited_haplotype/%s/haplotype_names.json" % (data["Name"], cell)
        d = json.load(open(jsonfile))
        if chrom not in d:
            continue
        names1 = d[chrom]["HP1"]
        names2 = d[chrom]["HP2"]
        print(cell, len(names1), len(names2), sep="\t")
        if len(names1) + len(names2) == 0:
            continue
        with pysam.AlignmentFile(bamfile) as f:
            if fw1 is None:
                header = f.header.as_dict()
                rgs = []
                for cell2 in cells:
                    rgs.append({'ID': cell2, 'SM': cell2, 'LB': cell2})
                header["RG"] = rgs
                fw1 = pysam.AlignmentFile(outfile1, "wb", header=header)
                fw2 = pysam.AlignmentFile(outfile2, "wb", header=header)
            for segment in f.fetch(chrom):
                dn = segment.get_tag("DN")
                fw = None
                if dn in names1:
                    fw = fw1
                elif dn in names2:
                    fw = fw2
                if fw is not None:
                    fw.write(segment)
    if fw1:
        fw1.close()
        fw2.close()

if __name__ == "__main__":
    main()
    