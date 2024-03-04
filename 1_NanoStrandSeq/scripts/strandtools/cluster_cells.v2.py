#!/usr/bin/env python
import sys, os, json
import numpy as np
import pandas as pd
import seaborn
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams["font.family"] = "arial"
import matplotlib.pyplot as plt
from PyPDF2 import PdfFileMerger
from pyBioInfo.IO.File import BedFile

# 专门为de novo组装的cluster_26改造的

def load_anchors(data):
    anchors = dict()
    for cell in data["Cells"]:
        path = data["Paths"][cell]
        # path = "results/%s/clusters/assigned/%s.bed.gz" % (data["Name"], cell)
        d = dict()
        for chrom in data["Chroms"]:
            d[chrom] = [dict(), dict()]
        with BedFile(path) as f:
            for obj in f:
                chrom = obj.chrom
                strand = obj.strand
                si = 0 if strand == "+" else 1
                start = obj.start
                base = obj.name[2]
                d[chrom][si][start] = base
        anchors[cell] = d
    return anchors
        
        
def compare(d1, d2):
    n1, n2 = 0, 0
    for p in d1.keys() & d2.keys():
        if d1[p] == d2[p]:
            n1 += 1
        else:
            n2 += 1
    return n1, n2
    
    
def main():
    infile, chroms, outdir = sys.argv[1:]
    chroms = chroms.split(",")
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    if not os.path.exists(outdir + "/meta"):
        os.mkdir(outdir + "/meta")
        
    # data = json.load(open(infile))
    data = dict()
    data["Cells"] = []
    data["Paths"] = dict()
    data["Chroms"] = chroms # ["cluster_26"]
    with open(infile) as f:
        for line in f:
            cell, path = line.strip("\n").split("\t")
            data["Cells"].append(cell)
            data["Paths"][cell] = path

    
    anchors = load_anchors(data)
    
    for chrom in data["Chroms"]:
        # if chrom != "cluster_26":
        #     continue
        select_cells = []
        for cell in data["Cells"]:
            d = anchors[cell][chrom]
            count = len(d[0]) + len(d[1])
            if count >= 300:
                select_cells.append(cell)
        select_cells.sort()
        print(chrom, len(select_cells), sep="\t")
        
        data1 = dict()
        data2 = dict()
        for c1 in select_cells:
            for c2 in select_cells:
                dc1, dw1 = anchors[c1][chrom]
                dc2, dw2 = anchors[c2][chrom]
                m = np.zeros((4, 2))
                m[0][0], m[0][1] = compare(dc1, dc2)
                m[1][0], m[1][1] = compare(dw1, dw2)
                m[2][0], m[2][1] = compare(dc1, dw2)
                m[3][0], m[3][1] = compare(dw1, dc2)
                m = [list(row) for row in m]
                data1["%s_%s" % (c1, c2)] = m
                n1 = m[0][0] + m[1][0] + m[2][1] + m[3][1] # CC-WW
                n2 = m[0][1] + m[1][1] + m[2][0] + m[3][0] # CW-WC
                data2[(c1, c2)] = [n1, n2]
        cells = select_cells
        with open(outdir + "/meta/%s.json" % chrom, "w+") as fw:
            json.dump(data1, fw)
        cluster1_cells = []
        cluster2_cells = []
        
        r = 1
        while len(cells) > 0:
            m = np.zeros((len(cells), len(cells)))
            for i, c1 in enumerate(cells):
                for j, c2 in enumerate(cells):
                    n1, n2 = data2[(c1, c2)]
                    log2fc = np.log2(np.divide(n1, n2))
                    v = 0
                    if max(n1, n2) >= 10:
                        v = log2fc
                        if v > 4:
                            v = 4
                        elif v < -4:
                            v = -4
                        elif np.isnan(v):
                            v = 0
                    m[i][j] = v
            d = pd.DataFrame(m)
            d.columns = cells
            d.index = cells

            ret = seaborn.clustermap(d, cmap="bwr", figsize=(8, 8), vmin=-1, vmax=1)
            ret.fig.savefig(outdir + "/meta/%s.R%d.clustermap.pdf" % (chrom, r), dpi=300)
            ret.data2d.to_csv(outdir + "/meta/%s.R%d.clustermap.tsv" % (chrom, r), sep="\t")
            plt.close()

            values = ret.data2d.values
            xs = []
            ys = []
            for n in range(0, len(values) + 1):
                s = 0
                for i in range(len(values)):
                    for j in range(len(values)):
                        v = int(values[i][j])
                        if i < n and j < n:
                            s += v
                        elif i >= n and j >= n:
                            s += v
                        else:
                            s -= v
                xs.append(n)
                ys.append(s)
            vmax = max(ys)
            plt.figure(figsize=(5, 3))
            plt.plot(xs, ys, marker="o")
            plt.tight_layout()
            plt.savefig(outdir + "/meta/%s.R%d.plot.pdf" % (chrom, r), dpi=300)
            plt.close()
            
            # pdf = PdfFileMerger()
            # for c in ret.data2d.columns:
            #     path = "results/%s/wc/%s/filtered.final/%s.pdf" % (data["Name"], c, chrom)
            #     pdf.append(path)
            # pdf.write(outdir + "/meta/%s.R%d.barplots.pdf" % (chrom, r))
                
            num_cluster1 = ys.index(max(ys))
            num_cluster2 = len(cells) - num_cluster1
            print("Number of cluster1:", num_cluster1)
            print("Number of cluster2:", num_cluster2)
            cluster1_cells = list(ret.data2d.columns[:num_cluster1])
            cluster2_cells = list(ret.data2d.columns[-num_cluster2:])

            rows = []
            for i, cell1 in enumerate(ret.data2d.columns):
                n1 = 0 # agree
                n2 = 0 # disagree
                for j, cell2 in enumerate(ret.data2d.columns):
                    v = ret.data2d.values[i][j]
                    if (i < num_cluster1) == (j < num_cluster1):
                        if v > 0:
                            n1 += 1
                        if v < 0:
                            n2 += 1
                    else:
                        if v > 0:
                            n2 += 1
                        if v < 0:
                            n1 += 1
                # print(n1, n2, cell1, sep="\t")
                rows.append([cell1, n1, n2])
            d2 = pd.DataFrame(rows)
            d2.columns = ["Cell", "Agree", "Disagree"]
            d2["Total"] = d2["Agree"] + d2["Disagree"]
            vmax = d2["Disagree"].max()
            if vmax / len(cells) < 0.2:
                break
            else:
                cells = list(sorted(d2[d2["Disagree"] < vmax]["Cell"]))
                r += 1
                continue
                
        with open(outdir + "/cells.%s.json" % chrom, "w+") as fw:
            d = {"Cluster1": cluster1_cells, "Cluster2": cluster2_cells}
            s = json.dumps(d, indent=4)
            fw.write(s)
            
    
if __name__ == '__main__':
    main()
    