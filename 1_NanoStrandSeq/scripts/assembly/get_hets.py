#!/usr/bin/env python
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
from pyBioInfo.IO.File import BedFile


def decode(s):
    ds = []
    if s == ".":
        return ds
    for items in s.split(";"):
        d = defaultdict(int)
        for item in items.split(","):
            k, v = item.split(":")
            v = int(v)
            d[k] = v
        ds.append(d)
    return ds


def get_score(base, ds):
    cell_count = 0
    score = 0
    if len(ds) >= 1:
        score = 0
        if len(ds) == 1:
            d = ds[0]
            v1 = d[base]
            v2 = sum(d.values())
            if v2 >= 2:
                cell_count += 1
            if v1 / v2 >= 0.8 and v2 >= 4:
                score = 1
        else:
            for d in ds:
                v1 = d[base]
                v2 = sum(d.values())
                if v2 >= 2:
                    cell_count += 1
                if v1 >= 2 and v1 / v2 >= 0.8:
                    score += 1
    return cell_count, score


def main():
    infile1, infile2, chrom, outfile = sys.argv[1:]
    
    hets = dict()
    with BedFile(infile1) as f:
        for obj in f:
            if obj.chrom == chrom:
                base1, base2 = obj.name.split("-")
                hets[obj.start] = [base1, base2]    
                
    matrix = np.zeros((4, 2), dtype=np.int)
    
    with open(infile2) as f, open(outfile, "w+") as fw:
        for line in f:
            row = line.strip("\n").split("\t")
            if row[7] == "." and row[8] == 0:
                continue
            pos = int(row[0])
            base_ref = row[1]
            base_pat, base_mat = row[2], row[3]
            hc = row[4] == "1"
            if base_pat == ".":
                base_pat = base_ref
            if base_mat == ".":
                base_mat = base_ref
            base_hp1, base_hp2 = row[5], row[6]
            ds1, ds2 = decode(row[7]), decode(row[8])
            cell_count1, score1 = get_score(base_hp1, ds1)
            cell_count2, score2 = get_score(base_hp2, ds2)
            conf1 = score1 >= 1 and score1 / cell_count1 >= 0.75
            conf2 = score2 >= 1 and score2 / cell_count2 >= 0.75
            # conf1 = score1 >= 2 and score1 / cell_count1 >= 0.75
            # conf2 = score2 >= 2 and score2 / cell_count2 >= 0.75
            
            keep = False
            if conf1:
                if conf2:
                    if base_hp1 != base_hp2:
                        keep = True
                else:
                    ret = hets.get(pos)
                    if ret is not None:
                        if base_hp1 == ret[0]:
                            base_hp2 = ret[1]
                            keep = True
                        elif base_hp1 == ret[1]:
                            base_hp2 = ret[0]
                            keep = True
            else:
                if conf2:
                    ret = hets.get(pos)
                    if ret is not None:
                        if base_hp2 == ret[0]:
                            base_hp1 = ret[1]
                            keep = True
                        elif base_hp2 == ret[1]:
                            base_hp1 = ret[0]
                            keep = True
                            
            if keep and base_hp1 != "-" and base_hp2 != "-":
                fw.write("\t".join(map(str, [chrom, pos, pos + 1, "%s-%s" % (base_hp1, base_hp2)])) + "\n")
                if hc:
                    if base_pat == base_hp1:
                        matrix[0][0] += 1
                    else:
                        matrix[0][1] += 1
                    if base_pat == base_hp2:
                        matrix[1][0] += 1
                    else:
                        matrix[1][1] += 1
                    if base_mat == base_hp1:
                        matrix[2][0] += 1
                    else:
                        matrix[2][1] += 1
                    if base_mat == base_hp2:
                        matrix[3][0] += 1
                    else:
                        matrix[3][1] += 1
                        
    d = pd.DataFrame(matrix)
    d.columns = ["Same", "Diff"]
    d.index = ["Pat_HP1", "Pat_HP2", "Mat_HP1", "Mat_HP2"]
    d["Precision"] = d["Same"] / d.sum(axis=1)
    d.index.name = "Comparison"
    d.to_csv(sys.stdout, sep="\t")
                            

if __name__ == '__main__':
    main()
    