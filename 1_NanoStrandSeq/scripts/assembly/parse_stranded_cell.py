#!/usr/bin/env python
import os
import sys
import re
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pyBioInfo.IO.File import BamFile


def get_edge(inliers):
    i1 = None
    i2 = None
    for i, flag in enumerate(inliers):
        if flag:
            if i1 is None:
                i1 = i
            i2 = i
    i2 += 1
    return i1, i2


def fit(points, results, min_num=50, cut=None):
    n = len(points)
    idxs = np.arange(len(points))
    if n < min_num:
        results.append(points)
    else:
        array = []
        xs = points[:, 0]
        ys = points[:, 1]
        for i in range(1000):
            idxs_selected = np.random.choice(idxs, 3, False)
            points_selected = points[idxs_selected]
            xs0 = points_selected[:, 0]
            ys0 = points_selected[:, 1]
            a, b = np.polyfit(xs0, ys0, 1)  
            v = np.sqrt(np.power(a, 2) + 1)
            ds = np.abs((ys - a * xs - b) / v)
            inliers = ds <= cut
            i1, i2 = get_edge(inliers)
            ds1 = ds[inliers]
            distance = sum(ds1)
            ad = np.mean(ds1)
            count = len(ds1)
            score = count ** 3 / distance
            cr = count / (i2 - i1)
            array.append([a, b, count, distance, ad, score, (i1, i2), cr])
        tmp = list(filter(lambda item: item[7] > 0.6, array))
        if len(tmp) > 0:
            array = list(sorted(tmp, key=lambda item: item[2]))
        else:
            array = list(sorted(array, key=lambda item: item[2]))
        item = array[-1]
        a, b = item[0], item[1]
        i1, i2 = item[6]
        results.append(points[i1:i2])
        if i1 > 0:
            assert len(points[0:i1])
            fit(points[0:i1], results, min_num, cut)
        if i2 < n:
            assert len(points[i2:])
            fit(points[i2:], results, min_num, cut)
           
            
def fits(xs, ys):
    np.random.seed(0)
    points = []
    for x, y in zip(xs, ys):
        points.append([x, y])
    points = np.array(points)
    results = list()
    min_point = max(50, int(len(points) * 0.02))
    threshold = max(20, int(len(points) * 0.002))
    fit(points, results, min_point, threshold)
    results = list(sorted(results, key=lambda item: item[0][0]))
    return results


def main():
    infile, outdir = sys.argv[1:]
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    alignments = []
    chroms = []
    sizes = []
    bam_header = None
    with BamFile(infile) as f:
        bam_header = f.handle.header.as_dict()
        for item in f.handle.header["SQ"]:
            chrom = item["SN"]
            size = item["LN"]
            if re.match("^chr([0-9]+|[XY])$", chrom):
                chroms.append(chrom)
                sizes.append(size) 
        for obj in f:
            alignments.append(obj)
        
    with open(outdir + "/reads.txt", "w+") as fw:
        fw.write("%d\n" % len(alignments))
        
    reads = defaultdict(list)
    for obj in alignments:
        reads[obj.chrom].append(obj)
        
    # dat1: CWR   
    rows = []
    for chrom, size in zip(chroms, sizes):
        t = 0
        c = 0
        w = 0
        c_p = 0
        c_m = 0
        w_p = 0
        w_m = 0
        for obj in reads[chrom]:
            try:
                xp = obj.segment.get_tag("XP")
            except Exception:
                xp = "U"
            t += 1
            if obj.strand == "+":
                c += 1
                if xp == "P":
                    c_p += 1
                elif xp == "M":
                    c_m += 1
            else:
                w += 1
                if xp == "P":
                    w_p += 1
                elif xp == "M":
                    w_m += 1
        rpm = t * 1e6 / size
        rows.append([chrom, size, t, rpm, c, c_p, c_m, w, w_p, w_m])
    dat1 = pd.DataFrame(rows)
    dat1.columns = ["Chrom", "Size", "Reads", "RPM", 
                    "Crick", "Crick.Paternal", "Crick.Maternal", 
                    "Watson", "Watson.Paternal", "Watson.Maternal"]
    dat1["CWR"] = np.log2(dat1["Crick"] / dat1["Watson"])
    dat1["CPS"] = np.log2((dat1["Crick.Paternal"] + dat1["Watson.Maternal"]) / (dat1["Crick.Maternal"] + dat1["Watson.Paternal"]))
    
    xs = np.abs(dat1["CWR"])
    ys = np.abs(dat1["CPS"])
    plt.figure(figsize=(4, 4))
    plt.scatter(xs, ys, marker=".")
    plt.xlabel("ABS(CWR)")
    plt.ylabel("ABS(CPS)")
    plt.tight_layout()
    plt.savefig(outdir + "/cwr_cps.raw.png", dpi=300)
    plt.close()
            
    # dat2: Bins
    WIDTH = 1000000
    rows = []
    for chrom, size in zip(chroms, sizes):
        regions = []
        for x in range(0, size, WIDTH):
            y = min(x + WIDTH, size)
            width = y - x
            if width < WIDTH / 2:
                regions[-1][1] = y
            else:
                regions.append([x, y])
        nbin = len(regions)
        widths = np.array([y - x for x, y in regions])
        covs = np.zeros(nbin) # all
        covs_c = np.zeros(nbin) # crick
        covs_w = np.zeros(nbin) # watson
        covs_c_p = np.zeros(nbin) # crick pat
        covs_c_m = np.zeros(nbin) # crick mat
        covs_w_p = np.zeros(nbin) # watson pat
        covs_w_m = np.zeros(nbin) # watson mat
        for obj in reads[chrom]:
            bi1 = int(obj.start / WIDTH)
            bi2 = int((obj.end - 1) / WIDTH)
            bi1 = min(bi1, nbin - 1)
            bi2 = min(bi2, nbin - 1)
            bi2 += 1
            xp = obj.segment.get_tag("XP")
            v = 1 / (bi2 - bi1)
            for bi in range(bi1, bi2):
                covs[bi] += v
                if obj.strand == "+":
                    covs_c[bi] += v
                    if xp == "P":
                        covs_c_p[bi] += v
                    elif xp == "M":
                        covs_c_m[bi] += v
                else:
                    covs_w[bi] += v
                    if xp == "P":
                        covs_w_p[bi] += v
                    elif xp == "M":
                        covs_w_m[bi] += v
        covs = covs * WIDTH / widths
        covs_c = covs_c * WIDTH / widths
        covs_w = covs_w * WIDTH / widths
        covs_c_p = covs_c_p * WIDTH / widths
        covs_c_m = covs_c_m * WIDTH / widths
        covs_w_p = covs_w_p * WIDTH / widths
        covs_w_m = covs_w_m * WIDTH / widths
        for i in range(nbin):
            region = regions[i]
            row = [chrom, i, region[0], region[1], widths[0], 
                covs[i], 
                covs_c[i], covs_c_p[i], covs_c_m[i], 
                covs_w[i], covs_w_p[i], covs_w_m[i]]
            rows.append(row)
    dat2 = pd.DataFrame(rows)
    dat2.columns = ["Chrom", "Bin", "Start", "End", "Width", "Total",
                    "Crick", "Crick.Paternal", "Crick.Maternal", 
                    "Watson", "Watson.Paternal", "Watson.Maternal"]
    dat2["CWR"] = np.log2(dat2["Crick"] / dat2["Watson"])

    nbin_max = max([len(dat2[dat2["Chrom"] == chrom]) for chrom in chroms])
    total_max = dat2["Total"].max()
    cw_max = max(dat2["Crick"].max(), dat2["Watson"].max())
    cov_mean = np.mean(dat2[dat2["Total"] > 0]["Total"])
    dat2["TooHigh"] = dat2["Total"] > cov_mean * 2
    cov_mean = np.mean(dat2[(dat2["Total"] > 0) & (~dat2["TooHigh"])]["Total"])     
    
    # scatter of cwr
    cs = []
    ws = []
    for chrom in chroms:
        d = dat2[(dat2["Chrom"] == chrom) & (~dat2["TooHigh"])]
        c = 0
        for v1, v2 in d[["Width", "Crick"]].values:
            c += v1 * v2 / WIDTH
        w = 0
        for v1, v2 in d[["Width", "Watson"]].values:
            w += v1 * v2 / WIDTH
        cs.append(c)
        ws.append(w)
    dat1["Crick.Filtered"] = cs
    dat1["Watson.Filtered"] = ws
    dat1["CWR.Filtered"] = np.log2(dat1["Crick.Filtered"] / dat1["Watson.Filtered"])
    dat1["CWR.Diff"] = dat1["CWR.Filtered"] - dat1["CWR"]
    vs1 = dat1["CWR"].abs()
    vs2 = dat1["CWR.Filtered"].abs()
    vmax = max(max(vs1), max(vs2))
    plt.figure(figsize=(5, 5))
    plt.scatter(vs1, vs2, marker="o", color="blue")
    plt.plot([-0.5, vmax + 0.5], [-0.5, vmax + 0.5], ls="--", color="grey")
    plt.xlabel("CWR of whole chromosome")
    plt.ylabel("CWR of whole chromosome (without too high bin)")
    plt.grid()
    plt.tight_layout()
    plt.savefig(outdir + "/cwr.filtered.png", dpi=300)
    
    dat1.to_csv(outdir + "/cwr.tsv", sep="\t", index=False)
    dat2.to_csv(outdir + "/bins.tsv", sep="\t", index=False)
    
    # barplot
    if not os.path.exists(outdir + "/barplot"):
        os.mkdir(outdir + "/barplot")
        
    for chrom, size in zip(chroms, sizes):
        d = dat2[dat2["Chrom"] == chrom]
        xs = d["Bin"].values
        
        fig, axs = plt.subplots(2, 2, figsize=(16, 4), sharex=True)
        plt.sca(axs[0][0])
        plt.title("%s (%s)" % (chrom, format(size, ",")))
        ys = d["Total"]
        colors = ["C0"] * len(xs)
        d1 = d[d["TooHigh"]]
        if len(d1):
            for x in d1["Bin"]:
                colors[x] = "red"
            
        plt.bar(xs, ys, width=1, color=colors)
        plt.plot([-0.5, xs[-1] + 0.5], [cov_mean, cov_mean], ls="--", lw=1, color="grey")
        plt.plot([-0.5, xs[-1] + 0.5], [0, 0], ls="-", lw=1, color="black", clip_on=False)
        plt.xlim(-0.5, nbin_max - 0.5)
        plt.ylim(0, total_max * 1.05)
        plt.ylabel("RPM")
        plt.tight_layout()
        
        plt.sca(axs[0][1])
        plt.title("%s (%s)" % (chrom, format(size, ",")))
        plt.bar(xs, ys, width=1, color=colors)
        plt.plot([-0.5, xs[-1] + 0.5], [cov_mean, cov_mean], ls="--", lw=1, color="grey")
        plt.plot([-0.5, xs[-1] + 0.5], [0, 0], ls="-", lw=1, color="black", clip_on=False)
        plt.xlim(-0.5, nbin_max - 0.5)
        plt.ylim(0, cov_mean * 2)
        plt.ylabel("RPM")
        plt.tight_layout()
        
        plt.sca(axs[1][0])
        ys1 = d["Crick"]
        ys2 = d["Crick.Paternal"]
        ys3 = d["Crick.Maternal"]
        ys4 = d["Watson"]
        ys5 = d["Watson.Paternal"]
        ys6 = d["Watson.Maternal"]
        plt.bar(xs, ys1, width=1, color="C0")
        plt.bar(xs, ys2, width=1, color="blue")
        plt.bar(xs, ys3, bottom=ys2, width=1, color="red")
        plt.bar(xs, -ys4, width=1, color="C1")
        plt.bar(xs, -ys6, width=1, color="red")
        plt.bar(xs, -ys5, bottom=-ys6, width=1, color="blue")
        plt.plot([-0.5, xs[-1] + 0.5], [0, 0], ls="-", lw=1, color="black", clip_on=False)
        plt.xlim(-0.5, nbin_max - 0.5)
        plt.ylim(-cw_max * 1.05, cw_max * 1.05)
        # plt.ylim(-cov_mean * 2, cov_mean * 2)
        plt.ylabel("RPM")
        plt.tight_layout()
        
        plt.sca(axs[1][1])
        plt.bar(xs, ys1, width=1, color="C0")
        plt.bar(xs, ys2, width=1, color="blue")
        plt.bar(xs, ys3, bottom=ys2, width=1, color="red")
        plt.bar(xs, -ys4, width=1, color="C1")
        plt.bar(xs, -ys6, width=1, color="red")
        plt.bar(xs, -ys5, bottom=-ys6, width=1, color="blue")
        plt.plot([-0.5, xs[-1] + 0.5], [0, 0], ls="-", lw=1, color="black", clip_on=False)
        plt.xlim(-0.5, nbin_max - 0.5)
        # plt.ylim(-cw_max * 1.05, cw_max * 1.05)
        plt.ylim(-cov_mean * 2, cov_mean * 2)
        plt.ylabel("RPM")
        plt.tight_layout()
        
    #     for ax in axs:
    #         ax.spines["top"].set_visible(False)
    #         ax.spines["right"].set_visible(False)
    #         ax.spines["bottom"].set_visible(False)
        
        plt.savefig(outdir + "/barplot/%s.png" % chrom, dpi=300)
        plt.close()
        
    # polyfit
    
    if not os.path.exists(outdir + "/fits"):
        os.mkdir(outdir + "/fits")
        
    array = []
    for chrom in chroms:
        ps = []
        ys = []
        y = 0
        for i, obj in enumerate(reads[chrom]):
            ps.append(obj.start)
            if obj.strand == "+":
                y += 1
            else:
                y -= 1
            ys.append(y)
        if len(ys) <= 100:
            array.append([])
            continue
        
        xs = np.arange(len(ys))
        ys = np.array(ys)
        ps = np.array(ps)
        ret = np.polyfit(xs, ys, 1)
        a, b = ret
        
        xs1 = np.array([xs[0], xs[-1]])
        ys1 = np.array([a * x + b for x in xs1])

        fig, axs = plt.subplots(1, 3, figsize=(12, 8))
        
        plt.sca(axs[0])
        plt.title(chrom)
        plt.plot(xs, ys, color="black")
        plt.xlim(xs[0], xs[-1])
        plt.ylim(-len(ys), len(ys))
        plt.tight_layout()
        
        plt.sca(axs[1])
        plt.title(chrom)
        plt.plot(xs, ys, color="blue")
        plt.plot(xs1, ys1, color="red", ls="--")
        plt.text(len(ys) * 0.1, len(ys) * 0.8, "%s" % np.poly1d(ret))
        plt.xlim(xs[0], xs[-1])
        plt.ylim(-len(ys), len(ys))
        plt.tight_layout()
        
        if len(ys) >= 100:
            results = fits(xs, ys)
        else:
            results = list()
        array.append(results)
        plt.sca(axs[2])
        plt.title(chrom)
        for i, ps in enumerate(results):
            xs1 = ps[:,0]
            ys1 = ps[:,1]
            a, b = 0, 0
            if len(xs1) >= 2:
                a, b = np.polyfit(xs1, ys1, 1)
            plt.plot(xs1, ys1, color="C%d" % i, lw=2, label="%.4f" % a)
        plt.xlim(xs[0], xs[-1])
        plt.ylim(-len(ys), len(ys))
        plt.legend()
        plt.tight_layout()
        
        plt.savefig(outdir + "/fits/%s.png" % chrom, dpi=300)
        plt.close()
        # break
        
    # max CWR
    part_chrom_conf_cwr = 0
    for chrom, size, results in zip(chroms, sizes, array):
        item = dat1[dat1["Chrom"] == chrom]
        rs = item["Reads"].values[0]
        rpm = item["RPM"].values[0]
        for i, points in enumerate(results):
            xs = points[:, 0]
            ys = points[:, 1]
            n = len(points)
            a, b = 0, 0
            if len(xs) >= 2:
                a, b = np.polyfit(xs, ys, 1)
            if rpm > 10 and n / rs >= 0.33 and abs(a) >= 0.9:
                cwr = abs(np.log2(np.divide(1 + a, 1 - a)))
                part_chrom_conf_cwr = max(part_chrom_conf_cwr, cwr)

    cwr_max = dat1["CWR.Filtered"].abs().max()
    if part_chrom_conf_cwr > 0:
        cwr_max = max(cwr_max, part_chrom_conf_cwr)
    with open(outdir + "/best_cwr.txt", "w+") as fw:
        fw.write("%f\n" % cwr_max)
        
    # output bam
    rows1 = []
    rows2 = []
    with BamFile(outdir + "/selected.bam", "wb", header=bam_header) as fw:
        for chrom, size, results in zip(chroms, sizes, array):
            item = dat1[dat1["Chrom"] == chrom]
            rs = item["Reads"].values[0]
            rpm = item["RPM"].values[0]
            selected = 0
            for i, points in enumerate(results):
                xs = points[:, 0]
                ys = points[:, 1]
                n = len(points)
                a, b = 0, 0
                if len(xs) >= 2:
                    a, b = np.polyfit(xs, ys, 1)
                if rpm > 10 and n / rs >= 0.33 and abs(a) >= 0.9:
                    cwr = abs(np.log2(np.divide(1 + a, 1 - a)))
                    part_chrom_conf_cwr = max(part_chrom_conf_cwr, cwr)
                wrote = False
                if cwr_max >= 4 and rpm > 10 and abs(a) < 0.3 and n / rs >= 0.02:
                    wrote = True
                    objs = reads[chrom]
                    for x in points[:, 0]:
                        obj = objs[x]
                        selected += 1
                        fw.write(obj)
                rows1.append([chrom, size, rs, rpm, i, points[0][0], points[-1][0], a, b, wrote])
            rows2.append([chrom, size, rs, rpm, selected])
    os.system("samtools sort -@ 4 -T %s -o %s %s" % (outdir + "/selected.sorted", outdir + "/selected.sorted.bam", outdir + "/selected.bam"))
    os.system("samtools index -@ 4 %s" % (outdir + "/selected.sorted.bam"))
    os.remove(outdir + "/selected.bam")
    
    if len(rows1) > 1:
        dat3 = pd.DataFrame(rows1)
        dat3.columns = ["Chrom", "Size", "Reads", "RPM", "LineID", "Index1", "Index2", "Slope", "Intercept", "Wrote"]
        dat3.to_csv(outdir + "/lines.tsv", sep="\t", index=False)
    if len(rows2) > 1:
        dat4 = pd.DataFrame(rows2)
        dat4.columns = ["Chrom", "Size", "Reads", "RPM", "Wrote"]
        dat4.to_csv(outdir + "/pass.tsv", sep="\t", index=False)
        
    print("Finished!")
    
            
if __name__ == '__main__':
    main()
    