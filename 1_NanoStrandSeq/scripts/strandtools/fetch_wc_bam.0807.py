#!/usr/bin/env python
import sys
import os
import numpy as np
import pysam
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt


def get_wc_duplicate_names(segments, outdir):
    names = []
    if len(segments) < 100:
        return names
    
    chrom = segments[0].reference_name
    
    x = 0
    y = 0
    xs = [0]
    ys = [0]
    for segment in segments:
        x += 1
        if segment.is_reverse:
            y -= 1
        else:
            y += 1
        xs.append(x)
        ys.append(y)

    w = int(len(xs) * 0.01)
    w = max(w, 10)
    ks = [] 
    for i1 in range(0, len(xs)):
        i2 = i1 + w
        if i2 >= len(xs):
            break
        v1, v2 = ys[i1], ys[i2]
        k = abs((v2 - v1) / w)
        ks.append(k)
    ks = np.array(ks)
    
    plt.figure(figsize=(4, 3))
    plt.plot(np.arange(len(ks)), ks)
    plt.ylim(-0.1, 1.1)
    plt.axhline(0.6, ls="--", lw=1, color="grey")
    plt.xlabel("Window index")
    plt.ylabel("|K|")
    plt.tight_layout()
    plt.savefig(outdir + "/windows.%s.png" % chrom, dpi=300)
    plt.close()
    
    # peaks
    peaks = [] 
    s = None
    for i, k in enumerate(ks):
        if k > 0.6:
            if s is None:
                s = i
        else:
            if s is not None:
                peaks.append([s, i])
                s = None
    if s is not None:
        peaks.append([s, i])
        
    # reads
    reads = [] 
    for s, e in peaks:
        reads.append([s, e + w])
        
    # extended reads
    tmp = []
    for s, e in reads:
        s = max(s - w, 0)
        e = min(e + w, len(segments))
        tmp.append([s, e])
    reads = tmp
    
    # merge read range
    i = 0
    while i < len(reads) - 1:
        if reads[i + 1][0] <= reads[i][1]:
            reads[i][1] = max(reads[i][1], reads[i + 1][1])
            reads.pop(i + 1)
        else:
            i += 1

    # exclude read range
    tmp = []
    if len(reads) == 0:
        tmp.append([0, len(segments)])
    else:
        e0 = 0
        for s, e in reads:
            if s > e0:
                tmp.append([e0, s])
            e0 = e
        if len(segments) > e0:
            tmp.append([e0, len(segments)])
    reads = tmp
            
    xlim1, xlim2 = 0, len(segments)
    ylim1 = (max(ys) + min(ys)) / 2 - (xlim2 - xlim1) / 2
    ylim2 = (max(ys) + min(ys)) / 2 + (xlim2 - xlim1) / 2
    plt.figure(figsize=(5, 5))
    plt.plot(xs, ys)
    for s, e in reads:
        plt.fill_betweenx([ylim1, ylim2], s, e, color="lightgrey")
    plt.xlim(xlim1, xlim2)
    plt.ylim(ylim1, ylim2)
    plt.xlabel("Read index")
    plt.tight_layout()
    plt.savefig(outdir + "/reads.%s.png" % chrom, dpi=300)
    plt.close()
    
    for s, e in reads:
        for segment in segments[s:e]:
            names.append(segment.get_tag("DN"))
    return names
    
def process_chrom(f, chrom, fw, outdir):
    segments1 = [] # all segments
    segments2 = [] # uniq segments
    for segment in f.fetch(chrom):
        segments1.append(segment)
        if segment.is_duplicate:
            continue
        if segment.get_tag("BH") == "Y":
            continue
        segments2.append(segment)
    names = get_wc_duplicate_names(segments2, outdir)
    n = 0
    for segment in segments1:
        if segment.get_tag("DN") in names:
            fw.write(segment)
            n += 1
    print(chrom, len(segments1), len(segments2), n, sep="\t")
    

def main():
    infile, outdir = sys.argv[1:]
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    print("Chrom\tReads\tUniqReads\tWCReads")
    with pysam.AlignmentFile(infile) as f, \
        pysam.AlignmentFile(outdir + "/wc.bam", "wb", f) as fw:
        for chrom in f.references:
            # if chrom == "chrY" or chrom == "chrM":
            #     continue
            process_chrom(f, chrom, fw, outdir)
            # break
    cmd = "samtools index %s" % (outdir + "/wc.bam")
    assert os.system(cmd) == 0
                    
    
if __name__ == '__main__':
    main()
    