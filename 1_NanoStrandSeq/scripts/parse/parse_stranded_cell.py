#!/usr/bin/env python
import sys
import os
import json
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import norm
import re
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from pyBioInfo.IO.File import BedFile, BamFile
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import ShiftLoader


def load_reads(data):
    chromosomes = dict()
    reads = defaultdict(list)
    chroms = []
    sam_header = None
    with BamFile(data["bam"]) as f:
        sam_header = f.handle.header.as_dict()
        for item in sam_header["SQ"]:
            chrom = item["SN"]
            length = item["LN"]
            if chrom == "chrM":
                continue
            chroms.append(chrom)
            chromosomes[chrom] = GRange(chrom=chrom, start=0, end=length, name=chrom)
        for obj in f:
            if obj.chrom == "chrM":
                continue
            obj.parental = obj.segment.get_tag("XP")
            reads[obj.chrom].append(obj)
    data["chroms"] = chroms
    data["chromosomes"] = chromosomes
    data["sam_header"] = sam_header
    data["reads"] = reads
    return data
        
        
def load_blanks(data):
    blanks = defaultdict(list)
    with BedFile(data["bed"]) as f:
        for obj in f:
            blanks[obj.chrom].append(obj)
    data["blanks"] = blanks
    return data


def mark_chromosome_blanks(data):
    for chrom in data["chroms"]:
        v = 0
        for b in data["blanks"][chrom]:
            v += len(b)
        c = data["chromosomes"][chrom]
        c.blank = v
        c.fill = len(c) - v
    return data


def report_chromosome_information(data):
    # table
    rows = []
    for chrom in data["chroms"]:
        c = data["chromosomes"][chrom]
        rows.append([c.name, len(c), c.blank, len(data["reads"][chrom])])
    dat = pd.DataFrame(rows)
    dat.columns = ["Chrom", "Length", "Blank", "Reads"]
    dat["BlankRatio"] = dat["Blank"] / dat["Length"]
    dat["Clean"] = dat["Length"] - dat["Blank"]
    dat["RPM"] = dat["Reads"] * 1e6 / dat["Clean"]
    dat["Autosomal"] = [re.match("^chr[0-9]+$", chrom) is not None for chrom in dat["Chrom"]]
    dat.to_csv("%s/chromosomes.info.tsv" % data["outdir"], sep="\t", index=False)
    
    # plot
    tmp = dat[dat["Autosomal"]]
    rpm_mean = tmp["RPM"].mean()
    rpm_std = tmp["RPM"].std()
    std_ratio = rpm_std / rpm_mean
    xs = np.arange(len(dat))
    ys = dat["RPM"]
    xticks = dat["Chrom"]
    plt.figure(figsize=(8, 3))
    plt.title("Cell: %s, Mean: %.2f, Std: %.2f, Ratio: %.4f" % (data["cell"], rpm_mean, rpm_std, std_ratio))
    colors = ["C0" if x else "C1" for x in dat["Autosomal"]]
    plt.bar(xs, ys, edgecolor="black", color=colors)
    plt.axhline(rpm_mean, color="C2")
    plt.xlim(min(xs) - 0.5, max(xs) + 0.5)
    plt.xticks(xs, xticks, rotation=45)
    plt.ylabel("Reads / Million Base Range")
    plt.tight_layout()
    plt.savefig("%s/chromosome.png" % data["outdir"], dpi=300)
    plt.close()
    
    # sex
    xr = dat[dat["Chrom"] == "chrX"]["RPM"].values[0] / rpm_mean
    yr = dat[dat["Chrom"] == "chrY"]["RPM"].values[0] / rpm_mean
    if 0.25 < xr < 0.75 and 0.25 < yr < 0.75:
        sex = "male"
    elif 0.75 < xr and yr < 0.25:
        sex = "female"
    else:
        sex = "unknown"
    data["chromosome_info"] = dat
    data["sex"] = sex
    return data


def make_bins(data):
    # Binning
    bin_width = data["bin_width"]
    bins = defaultdict(list)
    for chrom in data["chroms"]:
        rows = []
        length = len(data["chromosomes"][chrom])
        for start in np.arange(0, length, bin_width):
            end = min(start + bin_width, length)
            # The bin that is smaller than 0.25 * bin width 
            # will be merged with the previous bin.
            if end - start > bin_width * 0.25 or len(rows) == 0:
                rows.append([start, end])
            else:
                rows[-1][1] = end
        for i, (start, end) in enumerate(rows):
            obj = GRange(chrom=chrom, start=start, end=end)
            obj.idx = i
            bins[chrom].append(obj)
    data["bins"] = bins
    return data


def mark_bin_blanks(data):
    for chrom in data["chroms"]:
        loader = ShiftLoader(data["blanks"][chrom])
        for b in data["bins"][chrom]:
            v = 0
            for obj1 in loader.fetch(obj=b):
                start = max(b.start, obj1.start)
                end = min(b.end, obj1.end)
                assert start < end
                v += (end - start)
            b.blank = v
            b.fill = len(b) - v
    return data


def count_bin_reads(data):
    for chrom in data["chroms"]:
        loader = ShiftLoader(data["reads"][chrom])
        for b in data["bins"][chrom]:
            b.crick = 0
            b.watson = 0
            b.crick_paternal = 0
            b.crick_maternal = 0
            b.watson_paternal = 0
            b.watson_maternal = 0
            for align in loader.fetch(obj=b):
                if align.strand == "+":
                    b.crick += 1
                    if align.parental == "P":
                        b.crick_paternal += 1
                    elif align.parental == "M":
                        b.crick_maternal += 1
                else:
                    b.watson += 1
                    if align.parental == "P":
                        b.watson_paternal += 1
                    elif align.parental == "M":
                        b.watson_maternal += 1
            b.crick_norm = np.divide(b.crick * 1e6, b.fill) if b.fill > 0 else 0
            b.crick_paternal_norm = np.divide(b.crick_paternal * 1e6, b.fill) if b.fill > 0 else 0
            b.crick_maternal_norm = np.divide(b.crick_maternal * 1e6, b.fill) if b.fill > 0 else 0
            b.watson_norm = np.divide(b.watson * 1e6, b.fill) if b.fill > 0 else 0
            b.watson_paternal_norm = np.divide(b.watson_paternal * 1e6, b.fill) if b.fill > 0 else 0
            b.watson_maternal_norm = np.divide(b.watson_maternal * 1e6, b.fill) if b.fill > 0 else 0
    return data


def report_bin_information(data):
    rows = []
    for chrom in data["chroms"]:
        for b in data["bins"][chrom]:
            rows.append([b.chrom, b.start, b.end, len(b), b.idx, b.blank, b.fill, 
                         b.crick, b.crick_paternal, b.crick_maternal, 
                         b.watson, b.watson_paternal, b.watson_maternal, 
                         b.crick_norm, b.crick_paternal_norm, b.crick_maternal_norm, 
                         b.watson_norm, b.watson_paternal_norm, b.watson_maternal_norm])
    dat = pd.DataFrame(rows)
    dat.columns = ["Chrom", "Start", "End", "Length", "Index", "Blank", "Clean", 
                    "Crick", "Crick.Paternal", "Crick.Maternal", 
                    "Watson", "Watson.Paternal", "Watson.Maternal", 
                    "Crick.Norm", "Crick.Paternal.Norm", "Crick.Maternal.Norm", 
                    "Watson.Norm", "Watson.Paternal.Norm", "Watson.Maternal.Norm"]
    dat.to_csv("%s/bins.info.tsv" % data["outdir"], sep="\t", index=False)
    data["bin_info"] = dat
    return data


def plot_bins(data, outdir, show_parental=True, show_blanks=False, show_regions=False, show_partitions=False, show_outliers=True):
    # print(outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    sex = data["sex"]
    bin_width = data["bin_width"]
    chroms = data["chroms"]
    bins = data["bins"]
    blanks = data["blanks"]
    chromosomes = data["chromosomes"]
    reads = data["reads"]
    
    vs = []
    for chrom in chroms:
        if sex == "female" and chrom == "chrY":
            continue
        for b in bins[chrom]:
            if b.fill < 0.5 * bin_width:
                continue
            vs.append(b.crick_norm + b.watson_norm)
    bin_count_max = max([len(v) for v in bins.values()])
    mean = np.mean(vs)
    std = np.std(vs)
    ylim = mean * 3
    
    for chrom in chroms:
        bs = bins[chrom]
        xs = np.arange(len(bs))
        ys1 = np.array([b.crick_norm for b in bs])
        ys2 = np.array([b.crick_paternal_norm for b in bs])
        ys3 = np.array([b.crick_maternal_norm for b in bs])
        ys4 = np.array([b.watson_norm for b in bs])
        ys5 = np.array([b.watson_paternal_norm for b in bs])
        ys6 = np.array([b.watson_maternal_norm for b in bs])
    
        plt.figure(figsize=(12, 3))
        
        plt.bar(xs, ys1, width=1, color="C0", label="Crick", zorder=10)
        plt.bar(xs, -ys4, width=1, color="C1", label="Watson", zorder=10)
        if show_parental:
            plt.bar(xs, ys2, width=1, color="blue", label="Paternal", zorder=10)
            plt.bar(xs, ys3, bottom=ys2, width=1, color="red", label="Maternal", zorder=10)
            plt.bar(xs, -ys5, bottom=-ys6, width=1, color="blue", zorder=10)
            plt.bar(xs, -ys6, width=1, color="red", zorder=10)
        plt.xlim(0.5, bin_count_max - 0.5)
        plt.ylim(-ylim, ylim)
        plt.plot([-0.5, xs[-1] + 0.5], [0, 0], color="grey", lw=2, zorder=0)
        plt.ylabel(chrom)
        
        if show_outliers:
            xs0 = []
            ys0 = []
            for x, y1, y2 in zip(xs, ys1, ys4):
                if y1 >= ylim:
                    xs0.append(x)
                    ys0.append(ylim * 0.98)
                if y2 >= ylim:
                    xs0.append(x)
                    ys0.append(-ylim * 0.98)
            plt.scatter(xs0, ys0, marker="o", color="red", s=50, clip_on=False, zorder=20)
        
        if show_blanks:
            for start, end in blanks[chrom]:
                x1 = start / bin_width
                x2 = end / bin_width
                w = x2 - x1
                x = (x1 + x2) / 2
                plt.bar(x, ylim * 0.1, width=w, bottom=-ylim * 1, color="grey", lw=8)
                
        if show_regions:
            ci = 0
            for start, end in data["regions"][chrom]:
                x1 = start / bin_width
                x2 = end / bin_width
                w = x2 - x1
                x = (x1 + x2) / 2
                plt.bar(x, ylim * 0.1, width=w, bottom=ylim * 0.9, color="C%d" % (ci % 10), lw=8)
                ci += 1
                
        if show_partitions:
            ci = 0
            for obj in data["partitions"][chrom]:
                if abs(obj.region.line.a) > 0.3:
                    continue
                if obj.fill < chromosomes[chrom].fill * 0.1:
                    continue
                if len(obj.reads) < len(reads[chrom]) * 0.1:
                    continue
                x1 = obj.start / bin_width
                x2 = obj.end / bin_width
                w = x2 - x1
                x = (x1 + x2) / 2
                plt.bar(x, ylim * 0.1, width=w, bottom=ylim * 0.75, color="C%d" % (ci % 10), lw=8)
                ci += 1
        
                
        show_edge = False
        if not show_edge:
            ax = plt.gca()
            ax.spines["top"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.spines["right"].set_visible(False)

        plt.legend(ncol=1, loc="upper left", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        if outdir is None:
            plt.show()
        else:
            plt.savefig(outdir + "/%s.png" % chrom, dpi=300)
        plt.close()
        
class LineModel(object):
    def __init__(self, points, sampling, max_dis=None):
        self.a = 0
        self.b = 0
        self.width = None
        self.avg_dis = None
        self.points = None
        self.score = 0
        self.outliers = []
        xs = points[:, 0]
        ys = points[:, 1]
        if sampling and len(points) > 10:
            idxs = np.arange(len(points))
            idxs1 = np.random.choice(idxs, 3, False)
            points1 = points[idxs1]
            xs1 = points1[:, 0]
            ys1 = points1[:, 1]
            a1, b1 = np.polyfit(xs1, ys1, 1)
            v = np.sqrt(np.power(a1, 2) + 1)
            distances = np.abs((ys - a1 * xs - b1) / v)
            inliers = distances <= max_dis
            i1, i2 = None, None
            for i0, flag in enumerate(inliers):
                if flag:
                    if i1 is None:
                        i1 = i0
                    i2 = i0
            i2 = i2 + 1
            points2 = points[i1:i2]
            if i1 > 0:
                self.outliers.append(points[:i1])
            if i2 < len(points):
                self.outliers.append(points[i2:])
            self.width = i2 - i1
            xs2 = points2[:, 0]
            ys2 = points2[:, 1]
            a2, b2 = np.polyfit(xs2, ys2, 1)
            v2 = np.sqrt(np.power(a2, 2) + 1)
            distances2 = np.abs((ys2 - a2 * xs2 - b2) / v2)
            # print(i1, i2, distances2)
            self.avg_dis = np.mean(distances2)
            self.points = points2
            self.a = a2
            self.b = b2
        else:
            self.points = points
            if len(points) > 1:
                a, b = np.polyfit(points[:, 0], points[:, 1], 1)
            elif len(points) == 1:
                a = 0
                b = points[0][1]
            else:
                assert False
            v = np.sqrt(np.power(a, 2) + 1)
            distances = np.abs((ys - a * xs - b) / v)
            self.avg_dis = np.mean(distances)
            self.width = len(points)
            self.a, self.b = a, b
        k = 0.001
        self.score = k * self.width - self.avg_dis
            
        
    def __str__(self):
        s = "+" if self.b >=0 else "-"
        return "y = %.4f x %s %.4f" % (self.a, s, abs(self.b))
    
    @classmethod
    def merge_lines(cls, lines):
        points = np.concatenate([line.points for line in lines], axis=0)
        return LineModel(points, False)
        
    
class RANSAC(object):
    def __init__(self, points):
        xs, ys = np.array(points[0]), np.array(points[1])
        assert len(xs) == len(ys)
        self.points = points
        self.lines = list()
        self.number = 1000 # number of sampling
        self.min_points = max(50, int(len(points) * 0.02))
        self.max_dis = max(20, int(len(points) * 0.002))
        np.random.seed(0)
        if len(points) >= 10:
            self._fit(points)
        self.lines = list(sorted(self.lines, key=lambda item: item.points[0][0]))
        
    def _fit(self, points):
        if len(points) < self.min_points:
            self.lines.append(LineModel(points, False))
        else:
            lines = []
            for i in range(self.number):
                lines.append(LineModel(points, True, self.max_dis))
            best = None      

            if True:
                lines = list(sorted(lines, key=lambda item: item.score))
                best = lines[-1]
            else:
                for line in lines:
                    if line.avg_dis > self.max_dis * 0.5:
                        continue
                    if best is None:
                        best = line
                    else:
                        if line.width > best.width:
                            best = line
            assert best is not None
            self.lines.append(best)
            
            for points1 in best.outliers:
                self._fit(points1)
                

def fit_lines(data):
    lines = defaultdict(list)
    for chrom in data["chroms"]:
        if chrom != "chrY":
            xs = [0]
            ys = [0]
            y = 0
            for i, item in enumerate(data["reads"][chrom]):
                if item.strand == "+":
                    y += 1
                else:
                    y -= 1
                xs.append(i + 1)
                ys.append(y)
            points = np.array([xs, ys]).T
            ransac = RANSAC(points)
            for line in ransac.lines:
                lines[chrom].append(line)
    data["lines_fitted"] = lines
    return data


def plot_lines(data, key, outdir=None):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    for chrom in data["chroms"]:
        lines = data[key][chrom]
        if len(lines) == 0:
            continue
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
        plt.sca(axs[0])
        plt.title(chrom)
        for line in lines:
            plt.plot(line.points[:, 0], line.points[:, 1], label="%.4f" % line.a)
        for i, line in enumerate(lines):
            if i == 0:
                continue
            v = line.points[0][0]
            plt.axvline(v, ls="--", lw=1, color="grey")
        count = sum([len(line.points) for line in lines])
        y_min = min([min(line.points[:, 1]) for line in lines])
        y_max = max([max(line.points[:, 1]) for line in lines])
        y_med = (y_min + y_max) / 2
        ylim1, ylim2 = y_med - count / 2, y_med + count / 2
        plt.xlabel("Read index")
        plt.ylabel("Cumulative")
        plt.xlim(0, count)
        plt.ylim(ylim1, ylim2)
        plt.grid(axis="y", ls="--")
        plt.legend(loc="upper left", bbox_to_anchor=(1, 1), ncol=4)
        axs[1].set_visible(False)
        if outdir is None:
            plt.show()
        else:
            plt.savefig("%s/%s.png" % (outdir, chrom), dpi=300)
        plt.close()
        

def merge_short_lines(data):
    d = dict()
    for chrom in data["chroms"]:
        array = data["lines_fitted"][chrom].copy()
        while True:
            n = 0
            i = 0
            while i < len(array) - 1:
                line1 = array[i]
                line2 = array[i + 1]
                merge = None
                if line1.a >= 0.9 and line2.a >= 0.9:
                    line3 = LineModel.merge_lines([line1, line2])
                    if line3.a > 0.9:
                        merge = line3
                elif line1.a <= -0.9 and line2.a <= -0.9:
                    line3 = LineModel.merge_lines([line1, line2])
                    if line3.a < -0.9:
                        merge = line3
                elif (abs(line1.a) < 0.3 and abs(line2.a) < 0.5) or (abs(line1.a) < 0.5 and abs(line2.a) < 0.3):
                    line3 = LineModel.merge_lines([line1, line2])
                    if abs(line3.a) < 0.3:
                        merge = line3
                if merge is None:
                    i += 1
                else:
                    array[i] = merge
                    array.pop(i + 1)
                    n += 1
            if n == 0:
                break 
        d[chrom] = array
    data["lines_merged"] = d
    return data


def merge_final_lines(data):
    d = dict()
    for chrom in data["chroms"]:
        array = data["lines_merged"][chrom].copy()
        margin = len(data["reads"][chrom]) * 0.05
        while True:
            n = 0
            i = 0
            while i < len(array) - 1:
                line1 = array[i]
                if abs(line1.a) > 0.9 and line1.width > margin:
                    i += 1
                    continue
                line2 = array[i + 1]
                if abs(line2.a) > 0.9 and line2.width > margin:
                    i += 1
                    continue 
                merge = LineModel.merge_lines([line1, line2])
                array[i] = merge
                array.pop(i + 1)
                n += 1
            if n == 0:
                break   
        d[chrom] = array
    data["lines_final"] = d
    return data


def make_region_and_partition(data):
    regions = defaultdict(list)
    partitions = defaultdict(list)
    for chrom in data["chroms"]:
        array = data["lines_final"][chrom]
        reads = data["reads"][chrom]
        
        for i, line in enumerate(array):
            i1 = line.points[0][0]
            i2 = line.points[-1][0]
            if i2 - i1 < 10:
                continue
            read1 = reads[i1]
            read2 = reads[i2 - 1]
            start = min(read1.start, read1.start)
            end = max(read1.end, read2.end)
            region = GRange(chrom=chrom, start=start, end=end, name="Region")
            region.line = line
            region.reads = reads[i1:i2]
            regions[chrom].append(region)

            discard1 = int((end - start) * 0.01)
            discard2 = int((i2 - i1) * 0.01)
            start1 = start
            end1 = end
            if i != 0:
                s1 = start + discard1
                s2 = reads[i1 + discard2].start
                start1 = max(s1, s2)
            if i != len(array) - 1:
                e1 = end - discard1
                e2 = reads[i2 - discard2 - 1].end
                end1 = min(e1, e2)
            tmp = []
            for read in reads[i1:i2]:
                if start1 <= read.start and read.end <= end1:
                    tmp.append(read)
            if start1 < end1:
                p = GRange(chrom=chrom, start=start1, end=end1, name="Partition")
                p.region = region
                p.reads = tmp
                partitions[chrom].append(p)

    data["regions"] = regions
    data["partitions"] = partitions
    return data


def mark_region_and_parition_blanks(data):
    for chrom in data["chroms"]:
        loader = ShiftLoader(data["blanks"][chrom])
        for obj in data["regions"][chrom]:
            v = 0
            for obj1 in loader.fetch(obj=obj):
                s = max(obj.start, obj1.start)
                e = min(obj.end, obj1.end)
                assert s < e
                v += (e - s)
            obj.blank = v
            obj.fill = len(obj) - v

        loader = ShiftLoader(data["blanks"][chrom])
        for obj in data["partitions"][chrom]:
            v = 0
            for obj1 in loader.fetch(obj=obj):
                s = max(obj.start, obj1.start)
                e = min(obj.end, obj1.end)
                assert s < e
                v += (e - s)
            obj.blank = v
            obj.fill = len(obj) - v
    return data


def output_bam(data):
    sam_header = data["sam_header"]
    out_bam = "%s/partitions.bam" % data["outdir"]
    rows = []
    with BamFile(out_bam, "wb", header=data["sam_header"]) as fw:
        for chrom in data["chroms"]:        
            pi = 0
            for partition in data["partitions"][chrom]:
                region = partition.region
                line = region.line
                a = line.a
                if abs(a) > 0.2:
                    continue
                r1 = partition.fill / data["chromosomes"][chrom].fill
                if r1 < 0.1:
                    continue
                r2 = len(partition.reads) / len(data["reads"][chrom])
                if r2 < 0.05:
                    continue
                rows.append([chrom, partition.start, partition.end, pi, len(partition.reads), r1, r2]) 
                for read in partition.reads:
                    read.segment.set_tag("PI", pi)
                    fw.write(read)
                pi += 1
    cmd = "samtools index %s" % out_bam
    assert os.system(cmd) == 0
                
    dat = pd.DataFrame(rows)
    dat.columns = ["Chrom", "Start", "End", "Index", "Reads", "FillRatio", "ReadRatio"]
    dat.to_csv("%s/partitions.info.tsv" % data["outdir"], sep="\t", index=False)
    data["partition_info"] = dat
    return data


def get_max_cwr(data):
    vs = []
    for chrom in data["chroms"]:
        for region in data["regions"][chrom]:
            if region.fill < data["chromosomes"][chrom].fill * 0.5:
                continue
            if len(region.reads) < len(data["reads"][chrom]) * 0.5:
                continue
            vs.append(region.line.a)
    if len(vs) > 0:
        vs = np.array(vs)
        k = max(abs(vs))
        cwr = (1 + k) / (1 - k)
        print(cwr)
    else:
        cwr = 1
    data["log2cwr"] = np.log2(cwr)
    return data


def main():
    infile1, infile2, outdir = sys.argv[1:]
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    data = dict()
    data["bam"] = infile1
    data["bed"] = infile2
    data["outdir"] = outdir
    data["cell"] = infile1.split("/")[-1][:-4]    
    data["bin_width"] = 1000000
    
    load_reads(data)
    load_blanks(data)
    
    mark_chromosome_blanks(data)
    report_chromosome_information(data)
    
    make_bins(data)
    mark_bin_blanks(data)
    count_bin_reads(data)
    report_bin_information(data)
    plot_bins(data, outdir + "/bins", show_parental=False)
    plot_bins(data, outdir + "/bins.parental", show_parental=True)
    
    fit_lines(data)
    plot_lines(data, "lines_fitted", outdir + "/lines_fitted")
    merge_short_lines(data)
    plot_lines(data, "lines_merged", outdir + "/lines_merged")
    merge_final_lines(data)
    plot_lines(data, "lines_final", outdir + "/lines_final")
    
    make_region_and_partition(data)
    mark_region_and_parition_blanks(data)
    plot_bins(data, outdir + "/bins.detail", show_parental=True, show_blanks=True, show_regions=True, show_partitions=True)

    output_bam(data)
    get_max_cwr(data)
    
    infos = dict()
    infos["bam"] = data["bam"]
    infos["bed"] = data["bed"]
    infos["outdir"] = data["outdir"]
    infos["bin_width"] = data["bin_width"]
    infos["cell"] = data["cell"]
    infos["log2cwr"] = data["log2cwr"]
    infos["reads"] = sum([len(vs) for vs in data["reads"].values()])
    with open(outdir + "/infos.tsv", "w+") as fw:
        fw.write(json.dumps(infos, indent=4))
    
    
if __name__ == '__main__':
    main()
    