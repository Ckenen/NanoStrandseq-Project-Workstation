#!/usr/bin/env python
import sys
import os
import re
import numpy as np
import pandas as pd
import optparse
import logging
from collections import defaultdict
from scipy.stats import ttest_ind_from_stats
import pysam
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams["font.family"] = "arial"
import matplotlib.pyplot as plt
from PyPDF2 import PdfFileMerger


MIN_MAPPING_QUALITY = 30
IGNORE_BLACKHOLE =True
IGNORE_SECODARY = True
IGNORE_SUPPLEMENTARY = True
REFERENCE_PATTERN = "chr([0-9]+|[XY])"
PLOT_READ_COUNT_BARPLOT = False


def load_bam(path):
    header = None
    segments = []
    with pysam.AlignmentFile(path) as f:
        header = f.header.as_dict()
        for segment in f:
            if segment.is_unmapped:
                continue
            if segment.mapping_quality < MIN_MAPPING_QUALITY:
                continue
            if IGNORE_SECODARY and segment.is_secondary:
                continue
            if IGNORE_SUPPLEMENTARY and segment.is_supplementary:
                continue
            if IGNORE_BLACKHOLE and segment.get_tag("BH") == "Y":
                continue
            segments.append(segment)
    return header, segments


def make_directory(path):
    if os.path.exists(path):
        logging.info("directory %s exists" % path)
    else:
        logging.info("created directory: %s" % path)
        os.mkdir(path)


def make_matrix(segments, chrom_length, bin_width):
    bin_count = int(chrom_length / bin_width)
    if chrom_length % bin_width > 0:
        bin_count += 1
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
    return matrix


def make_cumulative_points(segments):
    x, y = 0, 0
    xs, ys = [x], [y]
    for segment in segments:
        x += 1
        if segment.is_reverse:
            y -= 1
        else:
            y += 1
        xs.append(x)
        ys.append(y)
    xs, ys = np.array(xs), np.array(ys)
    return xs, ys


def filter_reads_by_k(segments, chrom, window_size, cutoff, outdir):
    segments = segments.copy()
    r = 1 # round
    while len(segments) >= window_size:
        xs, ys = make_cumulative_points(segments)
        title = "%s, size: %d, k: %.2f, round: %d" % (chrom, window_size, cutoff, r)
        if True:
            xlim1, xlim2 = 0, len(segments)
            ylim1 = (max(ys) + min(ys)) / 2 - (xlim2 - xlim1) / 2
            ylim2 = (max(ys) + min(ys)) / 2 + (xlim2 - xlim1) / 2
            plt.figure(figsize=(4, 4))
            plt.title(title)
            plt.plot(xs, ys)
            plt.xlim(xlim1, xlim2)
            plt.ylim(ylim1, ylim2)
            plt.xlabel("Read index")
            plt.tight_layout()
            plt.savefig(outdir + "/%s.R%d.XY.pdf" % (chrom, r), dpi=300)
            plt.close()
        
        ks = [] 
        for i1 in range(0, len(xs)):
            i2 = i1 + window_size
            if i2 >= len(xs):
                break
            v1, v2 = ys[i1], ys[i2]
            k = abs((v2 - v1) / window_size)
            ks.append(k)
        ks = np.array(ks)
        
        if True:
            plt.figure(figsize=(4, 4))
            plt.title(title)
            plt.plot(np.arange(len(ks)), ks)
            plt.ylim(-0.1, 1.1)
            plt.axhline(cutoff, ls="--", lw=1, color="grey")
            plt.xlabel("Window index")
            plt.ylabel("|K|")
            plt.tight_layout()
            plt.savefig(outdir + "/%s.R%d.KS.pdf" % (chrom, r), dpi=300)
            plt.close()
        
        i1 = None
        tmp1 = []
        for i2, k in enumerate(ks):
            if k >= cutoff:
                if i1 is None:
                    i1 = i2
            else:
                if i1 is not None:
                    tmp1.append([i1, i2])
                    i1 = None
        if i1 is not None:
            tmp1.append([i1, len(ks)])
        if len(tmp1) == 0:
            break

        tmp2 = []
        for i1, i2 in tmp1:
            i2 += window_size
            if len(tmp2) == 0:
                tmp2.append([i1, i2])
            elif i1 <= tmp2[-1][1]:
                tmp2[-1][1] = max(tmp2[-1][1], i2)
            else:
                tmp2.append([i1, i2])

        e = 0
        tmp3 = []
        for i1, i2 in tmp2:
            if i1 > e:
                tmp3.append([e, i1])
            e = i2
        if len(segments) > e:
            tmp3.append([e, len(segments)])

        tmp4 = []
        for i1, i2 in tmp3:
            if i2 - i1 < window_size:
                continue
            tmp4.extend(segments[i1:i2])
        segments = tmp4

        r += 1
        
    return segments


def main():
    # Option parser
    parser = optparse.OptionParser(usage="%prog [options] input.bam outdir")
    parser.add_option("-w", "--bin-width", dest="bin_width", default=1000000)
    parser.add_option("-p", "--pattern", dest="reference_pattern")
    options, args = parser.parse_args()
    bin_width = options.bin_width
    if len(args) != 2:
        parser.print_help()
        exit(1)
    
    # Logging config        
    logging.basicConfig(level=logging.INFO, 
                        format="[%(asctime)s %(levelname)s] %(message)s")

    infile, outdir = args
    matrix_dir = outdir + "/matrix"
    filter_cc_dir = outdir + "/filtered.cc"
    filter_ccw_dir = outdir + "/filtered.ccw"
    filter_final_dir = outdir + "/filtered.final"
    logging.info("input bam: %s" % infile)
    logging.info("output directory: %s" % outdir)
    logging.info("bin width: %d" % bin_width)
    make_directory(outdir)
    make_directory(matrix_dir)
    make_directory(filter_cc_dir)
    make_directory(filter_ccw_dir)
    make_directory(filter_final_dir)
    
    logging.info("loading bam file")    
    header, segments = load_bam(infile)
    chroms = []
    lengths = dict()
    for item in header["SQ"]:
        chrom, length = item["SN"], item["LN"]
        if re.match(REFERENCE_PATTERN, chrom) is None:
            continue
        chroms.append(chrom)
        lengths[chrom] = length
    logging.info("loaded %d chromosomes, %d segments" % (len(chroms), len(segments)))
      
    chrom_segments = defaultdict(list)
    for segment in segments:
        chrom_segments[segment.reference_name].append(segment)
    
    chrom_matrixs = dict()
    chrom_segments_nodup = dict()
    for chrom in chroms:
        segments1 = chrom_segments[chrom]
        segments2 = list(filter(lambda item: not item.is_duplicate, segments1))
        chrom_segments_nodup[chrom] = segments2
        logging.info("chromosome: %s, segments: %d, not duplicate segments: %d" % (chrom, 
                                                                                   len(segments1), 
                                                                                   len(segments2)))
        matrix = make_matrix(segments2, lengths[chrom], bin_width)
        chrom_matrixs[chrom] = matrix
        dat = pd.DataFrame(matrix)
        dat.columns = ["Crick", "Crick.Paternal", "Crick.Maternal", 
                       "Watson", "Watson.Paternal", "Watson.Maternal"]
        dat.index.name = "#Bin"
        dat.to_csv(matrix_dir + "/%s.tsv" % chrom, sep="\t")
        
    xs = []
    xticks = []
    offset = 0
    read_counts = []
    seps = []
    colors = []
    for ci, chrom in enumerate(chroms):
        if ci != 0:
            seps.append(offset)
        matrix = chrom_matrixs[chrom]
        color = "C%d" % (ci % 10)
        for idx in np.arange(len(matrix)):
            read_count = matrix[idx][0] + matrix[idx][3]
            read_counts.append(read_count)
            colors.append(color)
        xs.append(offset + len(matrix) / 2)
        xticks.append(chrom)
        offset += len(matrix)
    read_counts = np.array(read_counts)
    total_bin_count = len(read_counts)
    logging.info("total bin count: %d" % total_bin_count)
    read_counts1 = read_counts[read_counts > 0]
    logging.info("%s bin with read count > 0" % len(read_counts1))
    bin_read_mean = np.mean(read_counts1)
    bin_read_median = np.median(read_counts1)
    bin_read_std = np.std(read_counts1, ddof=1)
    logging.info("bin read mean: %.2f" % bin_read_mean)
    logging.info("bin read median: %d" % bin_read_median)
    logging.info("bin read std: %.2f" % bin_read_std)
        
    if PLOT_READ_COUNT_BARPLOT:
        logging.info("plot read count barplot for all bins")
        plt.figure(figsize=((total_bin_count + 150) * 0.005, 3))
        plt.title("mean: %.2f, median: %d, std: %.2f" % (bin_read_mean, bin_read_median, bin_read_std))
        plt.bar(np.arange(len(read_counts)) + 0.5, read_counts, width=1, color=colors)
        for x in seps:
            plt.axvline(x, lw=1, ls="--", color="grey")
        plt.axhline(bin_read_median, ls="--", lw=1, color="red")
        plt.xticks(xs, xticks, rotation=45)
        plt.xlim(0, len(read_counts))
        plt.ylabel("Read count / bin")
        plt.tight_layout()
        plt.savefig(outdir + "/barplot.pdf", dpi=300)
        plt.close()
        
    window_size = max(int(bin_read_median * 2), 50)
    chrom_segments_nocc = dict() # not CC or WW
    chrom_matrixs_nocc = dict()
    for chrom in chroms:
        logging.info("removing CC or WW reads for %s" % chrom)
        ss = filter_reads_by_k(segments=chrom_segments_nodup[chrom], 
                               chrom=chrom, 
                               window_size=window_size, 
                               cutoff=0.75,
                               outdir=filter_cc_dir)
        chrom_segments_nocc[chrom] = ss
        chrom_matrixs_nocc[chrom] = make_matrix(ss, lengths[chrom], bin_width)
    
    window_size = max(int(bin_read_median * 30), 750)
    chrom_segments_noccw = dict() # not CC or WW
    chrom_matrixs_noccw = dict()
    for chrom in chroms:
        logging.info("removing CCW or CWW reads for %s" % chrom)
        ss = filter_reads_by_k(segments=chrom_segments_nocc[chrom], 
                                chrom=chrom, 
                                window_size=window_size, 
                                cutoff=0.25,
                                outdir=filter_ccw_dir)
        chrom_segments_noccw[chrom] = ss
        chrom_matrixs_noccw[chrom] = make_matrix(ss, lengths[chrom], bin_width)
        
    # filter extremely high read count bin
    chrom_segments_final = dict()
    chrom_matrixs_final = dict()
    for chrom in chroms:
        matrix1 = chrom_matrixs[chrom]
        vs = matrix1[:, 0] + matrix1[:, 3]
        flags = [1] * len(matrix1)
        for i1 in range(0, len(vs)):
            i2 = i1 + 10
            if i2 > len(vs):
                break
            mean1 = np.mean(vs[i1:i2])
            std1 = np.std(vs[i1:i2], ddof=1)
            fc =  mean1 / bin_read_mean
            p = ttest_ind_from_stats(bin_read_mean, bin_read_std, 10, mean1, std1, 10)[1]
            if fc > 2 and p < 0.01:
                for i3 in range(i1, i2):
                    flags[i3] = 0
        tmp = []
        for s in chrom_segments_noccw[chrom]:
            idx = int(s.reference_start / bin_width)
            if flags[idx] == 0:
                continue
            tmp.append(s)
        chrom_segments_final[chrom] = tmp
        chrom_matrixs_final[chrom] = make_matrix(tmp, lengths[chrom], bin_width)
        
    logging.info("plot crick/watson barplot")
    max_bin_read_count = 0
    max_bin_count = 0
    for chrom in chroms:
        matrix1 = chrom_matrixs[chrom]
        max_bin_read_count = max(max_bin_read_count, max(matrix1[:, 0]))
        max_bin_read_count = max(max_bin_read_count, max(matrix1[:, 3]))
        max_bin_count = max(max_bin_count, len(matrix1))
    
    for chrom in chroms:
        matrix1 = chrom_matrixs[chrom]
        matrix2 = chrom_matrixs_noccw[chrom]
        matrix3 = chrom_matrixs_final[chrom]
        bin_count = len(matrix1)
        xs = np.arange(bin_count)
        ylim = max(max(matrix1[:, 0]), max(matrix1[:, 3])) * 1.2
        
        fig, axs = plt.subplots(1, 3, figsize=(18, 2), sharex=True, sharey=True)
        
        plt.sca(axs[0])
        plt.title("All reads")
        plt.bar(xs, matrix1[:, 0], width=1, color="C0")
        plt.bar(xs, matrix1[:, 1], width=1, color="blue")
        plt.bar(xs, matrix1[:, 2], bottom=matrix1[:, 1], width=1, color="red")
        plt.bar(xs, -matrix1[:, 3], width=1, color="C1")
        plt.bar(xs, -matrix1[:, 4], bottom=-matrix1[:, 5], width=1, color="blue")
        plt.bar(xs, -matrix1[:, 5], width=1, color="red")
        plt.plot([0, bin_count], [0, 0], lw=0.5, color="black")
        plt.xlim(0, max_bin_count)
        plt.ylim(-ylim, ylim)
        plt.xlabel("Bins (%d)" % bin_width)
        plt.ylabel(chrom)
        plt.tight_layout()

        plt.sca(axs[1])
        plt.title("WC reads")
        plt.bar(xs, matrix1[:, 0], width=1, color="grey")
        plt.bar(xs, -matrix1[:, 3], width=1, color="grey")
        plt.bar(xs, matrix2[:, 0], width=1, color="C0")
        plt.bar(xs, matrix2[:, 1], width=1, color="blue")
        plt.bar(xs, matrix2[:, 2], bottom=matrix2[:, 1], width=1, color="red")
        plt.bar(xs, -matrix2[:, 3], width=1, color="C1")
        plt.bar(xs, -matrix2[:, 4], bottom=-matrix2[:, 5], width=1, color="blue")
        plt.bar(xs, -matrix2[:, 5], width=1, color="red")
        plt.plot([0, bin_count], [0, 0], lw=0.5, color="black")
        plt.xlim(0, max_bin_count)
        plt.ylim(-ylim, ylim)
        plt.xlabel("Bins (%d)" % bin_width)
        plt.tight_layout()
        
        plt.sca(axs[2])
        plt.title("Final reads")
        plt.bar(xs, matrix1[:, 0], width=1, color="grey")
        plt.bar(xs, -matrix1[:, 3], width=1, color="grey")
        plt.bar(xs, matrix3[:, 0], width=1, color="C0")
        plt.bar(xs, matrix3[:, 1], width=1, color="blue")
        plt.bar(xs, matrix3[:, 2], bottom=matrix3[:, 1], width=1, color="red")
        plt.bar(xs, -matrix3[:, 3], width=1, color="C1")
        plt.bar(xs, -matrix3[:, 4], bottom=-matrix3[:, 5], width=1, color="blue")
        plt.bar(xs, -matrix3[:, 5], width=1, color="red")
        plt.plot([0, bin_count], [0, 0], lw=0.5, color="black")
        plt.xlim(0, max_bin_count)
        plt.ylim(-ylim, ylim)
        plt.xlabel("Bins (%d)" % bin_width)
        plt.tight_layout()
        
        plt.savefig(filter_final_dir + "/%s.pdf" % chrom, dpi=300)
        plt.close()
        
    pdf = PdfFileMerger()
    for chrom in chroms:
        path = filter_final_dir + "/%s.pdf" % chrom
        pdf.append(path)
    pdf.write(outdir + "/all_chroms.pdf")

    logging.info("output CW reads to bam file")
    
    fw = pysam.AlignmentFile(outdir + "/wc.bam", "wb", header=header)
    fw1 = open(outdir + "/wc.tsv", "w+")
    fw1.write("Chrom\tReads\tUniqReads\tWCReads\tWroteReads\n")
    for chrom in chroms:
        segments1 = chrom_segments[chrom]
        segments2 = chrom_segments_nodup[chrom]
        segments3 = chrom_segments_final[chrom]
        names = []
        for s in segments3:
            names.append(s.get_tag("DN"))
        names = set(names)
        count = 0
        for s in segments1:
            if s.get_tag("DN") in names:
                fw.write(s)
                count += 1
        fw1.write("\t".join(map(str, [chrom, 
                                      len(segments1), 
                                      len(segments2), 
                                      len(segments3), 
                                      count])) + "\n")
    fw1.close()
    fw.close()
    cmd = "samtools index %s" % (outdir + "/wc.bam")
    assert os.system(cmd) == 0
    
    
if __name__ == '__main__':
    main()
    