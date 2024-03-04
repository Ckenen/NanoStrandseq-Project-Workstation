#!/usr/bin/env python
import sys
import json
import pandas as pd
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import ShiftLoader


def load_svs(path):
    svs = []
    with pysam.VariantFile(path) as f:
        for record in f:
            if record.contig == "chrY":
                continue
            svtype = record.info["SVTYPE"]
            if svtype != "DEL" and svtype != "INS":
                continue
#             if list(record.filter)[0] != "PASS":
#                 continue
            svlen = abs(record.info["SVLEN"])
            sv = GRange(chrom=record.contig, start=record.start, end=record.stop, name=record.id)
            sv.record = record
            sv.svtype = svtype
            sv.svlen = svlen
            svs.append(sv)
    return svs


def filter_sv_by_regions(svs, regions):
    svs1 = []
    loader = ShiftLoader(regions)
    for sv in svs:
        n = len(list(loader.fetch(obj=sv)))
        if n == 0:
            svs1.append(sv)
    return svs1


def get_recall(svs_ref, svs_que, svtype):
    svs_ref = list(filter(lambda sv: sv.svtype == svtype, svs_ref))
    svs_que = list(filter(lambda sv: sv.svtype == svtype, svs_que))
    n_hit = 0
    loader = ShiftLoader(svs_que)
    for sv in svs_ref:
        hit = False
        for sv2 in loader.fetch(chrom=sv.chrom, start=sv.start - 1000, end=sv.end + 1000):
            if sv.svtype == sv2.svtype and min(sv.svlen, sv2.svlen) >= max(sv.svlen, sv2.svlen) * 0.7:
                hit = True
                break
        if hit:
            n_hit += 1
    data = dict()
    data["Reference"] = len(svs_ref)
    data["Query"] = len(svs_que)
    data["Reference_Hit"] = n_hit
    data["Reference_Recall"] = data["Reference_Hit"] / data["Reference"]
    
    return data


def main():
    
    f_vcf1, f_vcf2, f_quant1, f_quant2, f_bed, min_query_cell, outfile = sys.argv[1:]
    
    min_query_cell = int(min_query_cell)

    # load SVs
    
    svs_ref = load_svs(f_vcf1)
    svs_que = load_svs(f_vcf2)
    
    # SV names
    
    max_length = 10000
    min_length = 50
    min_reads = 5
    min_freq = 0.2

    dat = pd.read_csv(f_quant1, sep="\t")
    dat = dat[dat["Length"] <= max_length]
    dat = dat[dat["Length"] >= min_length]
    dat = dat[dat["Chrom"] != "chrY"]
    dat = dat[dat["AgreeRead"] >= min_reads]
    dat = dat[(dat["AgreeRead"] / (dat["AgreeRead"] + dat["DisagreeRead"])) >= min_freq]
    dat_ref = dat
    names_ref = set(dat["Name"])

    dat = pd.read_csv(f_quant2, sep="\t")
    dat = dat[dat["Length"] <= max_length]
    dat = dat[dat["Length"] >= min_length]
    dat = dat[dat["Chrom"] != "chrY"]
    dat = dat[dat["AgreeRead"] >= min_reads]
    if min_query_cell > 1:
        dat = dat[dat["AgreeCell"] >= min_query_cell]
    dat = dat[(dat["AgreeRead"] / (dat["AgreeRead"] + dat["DisagreeRead"])) >= min_freq]
    dat_que = dat
    names_que = set(dat["Name"])
    
    # filter sv by names
    
    svs_ref_1 = list(filter(lambda sv: sv.name in names_ref, svs_ref))
    svs_que_1 = list(filter(lambda sv: sv.name in names_que, svs_que))
    
    # blacklist regions
    
    regions = []
    with open(f_bed) as f:
        for line in f:
            chrom, start, end = line.strip("\n").split("\t")
            start, end = int(start), int(end)
            regions.append(GRange(chrom=chrom, start=start, end=end))
    regions.sort()
    
    # filter sv by regions
    
    svs_ref_2 = filter_sv_by_regions(svs_ref_1, regions)
    svs_que_2 = filter_sv_by_regions(svs_que_1, regions)
    
    # benchmark
    
    ref = svs_ref_2
    que = svs_que_2
    
    data = dict()
    
    # deletion
    
    d1 = get_recall(ref, que, "DEL")
    d2 = get_recall(que, ref, "DEL")
    recall = d1["Reference_Recall"]
    precision = d2["Reference_Recall"]
    f1 = 2 * recall * precision / (recall + precision)
    data["Del_Recall"] = recall
    data["Del_Precision"] = precision
    data["Del_F1"] = f1
    data["Del_Detail"] = [d1, d2]

    # precision
    
    d1 = get_recall(ref, que, "INS")
    d2 = get_recall(que, ref, "INS")
    recall = d1["Reference_Recall"]
    precision = d2["Reference_Recall"]
    f1 = 2 * recall * precision / (recall + precision)
    data["Ins_Recall"] = recall
    data["Ins_Precision"] = precision
    data["Ins_F1"] = f1
    data["Ins_Detail"] = [d1, d2]
    
    # output
    
    with open(outfile, "w+") as fw:
        json.dump(data, fw, indent=4)
    
    
if __name__ == "__main__":
    main()