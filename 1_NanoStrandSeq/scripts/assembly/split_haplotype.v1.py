#!/usr/bin/env python
import sys
import os
import re
import logging
import gzip
import json
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import pysam
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams["font.family"] = "arial"
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from PyPDF2 import PdfFileMerger
from pyBioInfo.Range import GRange
from pyBioInfo.Utils import ShiftLoader, SegmentTools

def mkdir(path):
    if os.path.exists(path):
        logging.info("directory %s exists" % path)
    else:
        logging.info("created directory: %s" % path)
        os.mkdir(path)


def load_snvs(bedfile):
    chrom_snvs = defaultdict(list)
    with gzip.open(bedfile, "rt") as f:
        for line in f:
            chrom, start, end, name = line.strip("\n").split("\t")[:4]
            start, end = int(start), int(end)
            a1, a2 = name.split("|")
            if a1 == a2 or len(a1) > 1 or len(a2) > 1:
                continue
            snv = GRange(chrom=chrom, start=start, end=end)
            snv.allele1 = a1
            snv.allele2 = a2
            chrom_snvs[snv.chrom].append(snv)
    return chrom_snvs


def get_parentals(segment, snps):
    parentals = []
    if len(snps) > 0:
        parsed_cigar = SegmentTools.parse_cigar(segment)
        sequence = segment.query_sequence
        for item in snps:
            start = item.start
            base = None
            for block in parsed_cigar:
                if block[3][0] <= start < block[3][1]:
                    if block[0] == "M":
                        offset = start - block[3][0]
                        read_idx = block[2][0] + offset
                        base = sequence[read_idx]
                    elif block[0] == "D":
                        base = "-"
                    else:
                        assert False
                    break
            if base == "-":
                parentals.append("-")
            elif base == item.allele1:
                parentals.append("P")  # Paternal
            elif base == item.allele2:
                parentals.append("M")  # Maternal
            else:
                parentals.append("O")  # Other
    return parentals


def determine_parental(parentals):
    parental = "U"  # Unknown
    if len(parentals) > 0:
        items = list(sorted(Counter(parentals).items(),
                     key=lambda item: item[1]))
        if len(items) > 0:
            parental = items[-1][0]
            if (len(parental) > 1) and (items[-2][1] == items[-1][1]):
                parental = "A"  # Ambigous
        else:
            assert False
    return parental


def mark_haplotype(chroms, chrom_segments, chrom_snvs):
    for chrom in chroms:
        segments = chrom_segments[chrom]
        snvs = chrom_snvs[chrom]
        loader = ShiftLoader(snvs) # SNV loader
        for segment in segments:
            start = segment.reference_start
            end = segment.reference_end
            snps = list(loader.fetch(chrom=chrom, start=start, end=end))
            parentals = get_parentals(segment, snps)
            parental = determine_parental(parentals)
            segment.set_tag("XP", parental)


def cal_bin_reads(chrom, length, segments, bin_width):
    bin_count = int(length / bin_width)
    if length % bin_width > 0:
        bin_count += 1
    rows = []
    for bi in range(bin_count):
        start = bi * bin_width
        end = min((bi + 1) * bin_width, length)
        rows.append([chrom, bi, start, end])
    d1 = pd.DataFrame(rows)
    d1.columns = ["Chrom", "Bin", "Start", "End"]
    matrix = np.zeros((bin_count, 6), dtype=np.int)
    for segment in segments:
        idx = int(segment.reference_start / bin_width)
        parental = segment.get_tag("XP")
        if segment.is_reverse: # -, watson
            matrix[idx][3] += 1
            if parental == "P":
                matrix[idx][4] += 1
            elif parental == "M":
                matrix[idx][5] += 1
        else:
            matrix[idx][0] += 1
            if parental == "P":
                matrix[idx][1] += 1
            elif parental == "M":
                matrix[idx][2] += 1
    d2 = pd.DataFrame(matrix)
    d2.columns = ["Crick", "Crick.P", "Crick.M", "Watson", "Watson.P", "Watson.M"]
    d2["Parental"] = d2[["Crick.P", "Crick.M", "Watson.P", "Watson.M"]].sum(axis=1)
    d = pd.concat([d1, d2], axis=1)
    return d

   
def main():
    bamfile, hetfile, outdir = sys.argv[1:]

    chrom_hp_duplicate_set_names = dict()
    
    # Logging config        
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s %(levelname)s] %(message)s")
    
    cell = os.path.splitext(os.path.basename(bamfile))[0]
    count_dir = outdir + "/counts"
    pdf_dir = outdir + "/pdfs"
    mkdir(outdir)
    mkdir(count_dir)
    mkdir(pdf_dir)

    # load bam
    
    logging.info("loading %s" % bamfile)

    chroms = []
    chrom_lengths = dict()
    chrom_segments = dict()
    with pysam.AlignmentFile(bamfile) as f:
        for chrom in f.references:
            if re.match("^chr([0-9]+|[X])$", chrom) is None:
                continue
            length = f.get_reference_length(chrom)
            chroms.append(chrom)
            chrom_lengths[chrom] = length
        for chrom in chroms:
            chrom_segments[chrom] = [s for s in f.fetch(chrom)]

    chrom_segments_rmdup = dict()
    for chrom, segments in chrom_segments.items():
        chrom_segments_rmdup[chrom] = list(filter(lambda s: not s.is_duplicate, segments))

    logging.info("loaded %d chromosomes" % len(chroms))
    logging.info("loaded %d segments" % sum([len(ss) for ss in chrom_segments.values()]))
    logging.info("-" * 80)
    logging.info("Chrom\tLength\tSegments\tRmDup")
    logging.info("-" * 80)
    for chrom in chroms:
        length = chrom_lengths[chrom]
        count1 = len(chrom_segments[chrom])
        count2 = len(chrom_segments_rmdup[chrom])
        logging.info("%s\t%d\t%d\t%d" % (chrom, length, count1, count2))
    logging.info("-" * 80)

    # mark haplotype
    chrom_snvs = load_snvs(hetfile)
    mark_haplotype(chroms, chrom_segments, chrom_snvs)
        
    # bin read counts
    
    logging.info("Calculating bin read counts")
    bin_width = 1000000
    logging.info("bin width is %d" % bin_width)
    chrom_counts = dict()
    chrom_counts_rmdup = dict()
    for chrom in chroms:
        logging.info("processing %s" % chrom)
        counts1 = cal_bin_reads(chrom, chrom_lengths[chrom], chrom_segments[chrom], bin_width)
        chrom_counts[chrom] = counts1
        counts1.to_csv(count_dir + "/%s.tsv" % chrom, sep="\t", index=False)
        counts2 = cal_bin_reads(chrom, chrom_lengths[chrom], chrom_segments_rmdup[chrom], bin_width)
        chrom_counts_rmdup[chrom] = counts2
        counts2.to_csv(count_dir + "/%s.rmdup.tsv" % chrom, sep="\t", index=False)
    d1 = pd.concat([chrom_counts[chrom] for chrom in chroms], axis=0, ignore_index=True)
    d2 = pd.concat([chrom_counts_rmdup[chrom] for chrom in chroms], axis=0, ignore_index=True)
    d1.to_csv(count_dir + "/all_chroms.tsv", sep="\t", index=False)
    d2.to_csv(count_dir + "/all_chroms.rmdup.tsv", sep="\t", index=False)
    
    # split haplotype
    
    max_bin_count = max([len(chrom_counts[chrom]) for chrom in chroms])
    logging.info("max bin count: %d" % max_bin_count)
    
    vs = d2["Parental"].values
    vs = vs[vs > 0]
    median = np.median(vs)
    window_size = max(int(median * 10), 10)
    logging.info("window size is %d" % window_size)
    
    logging.info("splitting haplotype")
    
    data = [] # [chrom, segments of hp1, segments of hp2]
    for chrom in chroms:
        # if chrom != "chr11":
        #     continue
        
        logging.info("-" * 80)
        logging.info("processing %s" % chrom)
        segments_crick = []
        segments_watson = []
        for segment in chrom_segments_rmdup[chrom]:
            if segment.is_reverse:
                segments_watson.append(segment)
            else:
                segments_crick.append(segment)
        logging.info("chrom: %s, crick segments: %d, watson segments: %d" % (chrom, len(segments_crick), len(segments_watson)))
        
        segments_selected = []
        segments_hp1 = []
        segments_hp2 = []
        
        fig = plt.figure(figsize=(14, 6))
        gs = gridspec.GridSpec(6, 14)
        plt.suptitle("%s.%s" % (cell, chrom))
        
        for strand_idx in range(2):
            if strand_idx == 0:
                segments_strand = segments_crick
                strand_name = "crick"
                ax1 = fig.add_subplot(gs[0:3, 0:3])
                ax2 = fig.add_subplot(gs[0:3, 3:6])
            else:
                segments_strand = segments_watson
                strand_name = "watson"
                ax1 = fig.add_subplot(gs[3:6, 0:3])
                ax2 = fig.add_subplot(gs[3:6, 3:6])
            
            x, y = 0, 0
            xs, ys = [x], [y]
            segments_parental = []
            for segment in segments_strand:
                parental = segment.get_tag("XP")
                if parental == "P":
                    x += 1
                    y += 1
                elif parental == "M":
                    x += 1
                    y -= 1
                else:
                    continue
                segments_parental.append(segment)
                xs.append(x)
                ys.append(y)
            logging.info("loaded %d segments with parental infomation for %s" % (len(segments_parental), strand_name))
            logging.info("generated %d points" % len(xs))
            
            ymed = (max(ys) + min(ys)) / 2
            ydiff = len(ys) / 2
            ylim1, ylim2 = ymed - ydiff, ymed + ydiff
            
            ks = []
            for i1 in range(len(ys) - window_size):
                i2 = i1 + window_size
                y1, y2 = ys[i1], ys[i2]
                k = abs((y2 - y1) / window_size)
                ks.append(k)
            logging.info("generated %d Ks" % len(ks))
                
            i1 = None
            array1 = [] # regions of K >= 0.8
            for i2, k in enumerate(ks):
                if k >= 0.8:
                    if i1 is None:
                        i1 = i2
                else:
                    if i1 is not None:
                        array1.append([i1, i2])
                        i1 = None
            if i1 is not None:
                array1.append([i1, len(ks)])
            logging.info("retain ks: %s" % str(array1))
            
            array2 = [] # regions of points
            for x1, x2 in array1:
                if x1 != 0:
                    x1 += int(window_size * 0.1)
                if x2 < len(ks):
                    x2 -= int(window_size * 0.1)
                x2 += window_size
                if x2 - x1 >= window_size * 0.5:
                    if len(array2) > 0 and x1 <= array2[-1][1]:
                        array2[-1][1] = max(array2[-1][1], x2)
                    else:
                        array2.append([x1, x2])
            logging.info("retain points: %s" % str(array2))
            
            for x1, x2 in array2:
                x2 -= 1
                y1, y2 = ys[x1], ys[x2]
                k = (y2 - y1) / (x2 - x1)
                logging.info("line: %d-%d, k: %.2f" % (x1, x2, k))
                segments_hp = None
                if k >= 0.8:
                    segments_hp = segments_hp1
                elif k <= -0.8:
                    segments_hp = segments_hp2
                if segments_hp is None:
                    continue
                assert x1 >= 0
                assert x2 <= len(segments_parental)
                segments_parental_selected = segments_parental[x1:x2]
                if x1 == 0:
                    start_min = 0
                else:
                    start_min = segments_parental_selected[0].reference_start
                if x2 == len(segments_parental):
                    start_max = chrom_lengths[chrom]
                else:
                    start_max = segments_parental_selected[-1].reference_start
                logging.info("selected %d parental segments, min start: %d, max start: %d" % (len(segments_parental_selected), start_min, start_max))
                segments_strand_selected = []
                for segment in segments_strand:
                    if segment.reference_start >= start_min and segment.reference_end <= start_max:
                        segments_strand_selected.append(segment)
                        segments_hp.append(segment)
                        segments_selected.append(segment)
                logging.info("selected %d segments" % len(segments_strand_selected))
            
            plt.sca(ax1)
            plt.title(strand_name)
            plt.plot(xs, ys)
            for x1, x2 in array2:
                plt.axvspan(x1, x2, color="lightgrey")
            if len(xs) > 0:
                plt.xlim(0, max(xs))
            if ylim1 < ylim2:
                plt.ylim(ylim1, ylim2)
            plt.xlabel("Read index (parental)")
            plt.ylabel("Cumulative number")
            plt.grid(ls="--")
            
            plt.sca(ax2)
            plt.title(strand_name)
            plt.plot(np.arange(len(ks)), ks)
            if len(ks) > 0:
                plt.xlim(0, len(ks))
            plt.ylim(-0.1, 1.1)
            for x1, x2 in array1:
                plt.axvspan(x1, x2, color="lightgrey")
            plt.xlabel("Window index")
            plt.ylabel("|K|")
            plt.grid(ls="--", color="grey")
        
        logging.info("plotting barplot")
        
        # print(len(segments_selected), len(segments_hp1), len(segments_hp2))
        d1 = chrom_counts_rmdup[chrom]
        d2 = cal_bin_reads(chrom, chrom_lengths[chrom], segments_selected, bin_width)
        d3 = cal_bin_reads(chrom, chrom_lengths[chrom], segments_hp1, bin_width)
        d4 = cal_bin_reads(chrom, chrom_lengths[chrom], segments_hp2, bin_width)
        ys1 = np.arange(len(d1)) + 0.5
        vs = d1["Crick"] + d1["Watson"]
        xlim = np.mean(vs) * 1.5
        # xlim = max(d1["Crick"].max(), d1["Watson"].max(), 10) * 1.2
        
        axs = [
            fig.add_subplot(gs[0:6, 6:8]),
            fig.add_subplot(gs[0:6, 8:10]),
            fig.add_subplot(gs[0:6, 10:12]),
            fig.add_subplot(gs[0:6, 12:14]),
        ]
        
        plt.sca(axs[0])
        plt.title("All")
        plt.barh(ys1, d1["Crick"], height=1, color="C0")
        plt.barh(ys1, d1["Crick.P"], height=1, color="blue")
        plt.barh(ys1, d1["Crick.M"], left=d1["Crick.P"], height=1, color="red")
        plt.barh(ys1, -d1["Watson"], height=1, color="C1")
        plt.barh(ys1, -d1["Watson.P"], left=-d1["Watson.M"], height=1, color="blue")
        plt.barh(ys1, -d1["Watson.M"], height=1, color="red")
        plt.ylabel("Bins")
        
        plt.sca(axs[1])
        plt.title("Selected")
        plt.barh(ys1, d1["Crick"], height=1, color="grey")
        plt.barh(ys1, -d1["Watson"], height=1, color="grey")
        plt.barh(ys1, d2["Crick"], height=1, color="C0")
        plt.barh(ys1, d2["Crick.P"], height=1, color="blue")
        plt.barh(ys1, d2["Crick.M"], left=d2["Crick.P"], height=1, color="red")
        plt.barh(ys1, -d2["Watson"], height=1, color="C1")
        plt.barh(ys1, -d2["Watson.P"], left=-d2["Watson.M"], height=1, color="blue")
        plt.barh(ys1, -d2["Watson.M"], height=1, color="red")
        
        plt.sca(axs[2])
        plt.title("Haplotype 1")
        plt.barh(ys1, d1["Crick"], height=1, color="grey")
        plt.barh(ys1, -d1["Watson"], height=1, color="grey")
        plt.barh(ys1, d3["Crick"], height=1, color="C0")
        plt.barh(ys1, d3["Crick.P"], height=1, color="blue")
        plt.barh(ys1, d3["Crick.M"], left=d3["Crick.P"], height=1, color="red")
        plt.barh(ys1, -d3["Watson"], height=1, color="C1")
        plt.barh(ys1, -d3["Watson.P"], left=-d3["Watson.M"], height=1, color="blue")
        plt.barh(ys1, -d3["Watson.M"], height=1, color="red")
        
        plt.sca(axs[3])
        plt.title("Haplotype 2")
        plt.barh(ys1, d1["Crick"], height=1, color="grey")
        plt.barh(ys1, -d1["Watson"], height=1, color="grey")
        plt.barh(ys1, d4["Crick"], height=1, color="C0")
        plt.barh(ys1, d4["Crick.P"], height=1, color="blue")
        plt.barh(ys1, d4["Crick.M"], left=d4["Crick.P"], height=1, color="red")
        plt.barh(ys1, -d4["Watson"], height=1, color="C1")
        plt.barh(ys1, -d4["Watson.P"], left=-d4["Watson.M"], height=1, color="blue")
        plt.barh(ys1, -d4["Watson.M"], height=1, color="red")
        
        for ax in axs:
            ax.set_xlim(-xlim, xlim)
            ax.set_ylim(-1, max_bin_count + 1)
            ax.plot([0, 0], [0, len(d1)], lw=0.5, color="black", zorder=999)
            # ax.spines["top"].set_visible(False)
            # ax.spines["left"].set_visible(False)
            # ax.spines["right"].set_visible(False)
        plt.tight_layout(rect=(0, 0, 1, 0.96))
        plt.savefig(pdf_dir + "/%s.pdf" % chrom, dpi=300)
        plt.close()
        
        # report
        logging.info("output to bam file")
        names1 = [] # duplicate set name
        names2 = []
        for segment in segments_hp1:
            names1.append(segment.get_tag("DN"))
        for segment in segments_hp2:
            names2.append(segment.get_tag("DN"))
        names1 = set(names1)
        names2 = set(names2)
        items1 = []
        items2 = []
        for segment in chrom_segments[chrom]:
            dn = segment.get_tag("DN")
            if dn in names1:
                items1.append(segment)
            elif dn in names2:
                items2.append(segment)
        data.append([chrom, len(items1), len(items2)])
        chrom_hp_duplicate_set_names[chrom] = dict()
        chrom_hp_duplicate_set_names[chrom]["HP1"] = list(names1)
        chrom_hp_duplicate_set_names[chrom]["HP2"] = list(names2)

    with open(outdir + "/haplotype_names.json", "w+") as fw:
        json.dump(chrom_hp_duplicate_set_names, fw)
    
    dat = pd.DataFrame(data)
    dat.columns = ["Chrom", "HP1", "HP2"]
    dat.to_csv(outdir + "/stats.tsv", sep="\t", index=False)
      
    logging.info("mergeing pdfs")
    pdf = PdfFileMerger()
    for chrom in chroms:
        path = pdf_dir + "/%s.pdf" % chrom
        pdf.append(path)
    pdf.write(outdir + "/all_chroms.pdf")
    
    logging.info("completed!")
    

if __name__ == '__main__':
    main()
    