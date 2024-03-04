#!/usr/bin/env python
import sys
import json
import pysam
        
        
def main():
    f_cluster, assembly_dir, chrom, outfile1, outfile2 = sys.argv[1:]
    
    with open(f_cluster) as f:
            d = json.load(f)
            cluster1 = d["Cluster1"]
            cluster2 = d["Cluster2"]
            
    fw1, fw2 = None, None
    cells2 = cluster1 + cluster2

    for cell in cells2:
        bamfile = assembly_dir + "/prepare/bams/%s.bam" % cell
        jsonfile = assembly_dir + "/wc/%s/duplicate_set_names.json" % cell
        names = set(json.load(open(jsonfile))[chrom]) # dupliate set names
        #print(cell, len(names), sep="\t")
        if len(names) == 0:
            continue
        with pysam.AlignmentFile(bamfile) as f:
            if fw1 is None:
                header = f.header.as_dict()
                rgs = []
                for cell2 in cells2:
                    rgs.append({'ID': cell2, 'SM': cell2, 'LB': cell2})
                header["RG"] = rgs
                fw1 = pysam.AlignmentFile(outfile1, "wb", header=header)
                fw2 = pysam.AlignmentFile(outfile2, "wb", header=header)
            w1, w2 = fw1, fw2
            if cell in cluster2:
                w1, w2 = fw2, fw1
            for segment in f.fetch(chrom):
                dn = segment.get_tag("DN")
                if dn not in names:
                    continue
                if segment.is_reverse:
                    w = w2
                else:
                    w = w1
                w.write(segment)
    if fw1:
        fw1.close()
        fw2.close()

if __name__ == "__main__":
    main()
    
