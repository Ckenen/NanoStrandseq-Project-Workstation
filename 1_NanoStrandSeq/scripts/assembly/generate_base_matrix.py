#!/usr/bin/env python
import sys
import os
import pysam
import multiprocessing
import gzip
from collections import defaultdict, Counter
import numpy as np
from pyBioInfo.IO.File import BamFile, Alignment
from pyBioInfo.Utils import BundleBuilder, SegmentTools


def load_alignments(path, chrom, start, end):
    with pysam.AlignmentFile(path) as f:
        for segment in f.fetch(chrom, start, end):
            obj = Alignment(segment)
            events = SegmentTools.get_events(segment)
            obj.me = events
            yield obj
            
            
def execute_bin(task):
    infile1 = task["bam"]
    infile2 = task["fasta"]
    chrom = task["chrom"]
    bin_start = task["start"]
    bin_end = task["end"]
    outfile = task["outfile"]
    
    sequence = ""
    with pysam.FastaFile(infile2) as fasta:
        sequence = fasta.fetch(chrom, bin_start, bin_end)
    
    if outfile.endswith(".gz"):
        fw = gzip.open(outfile, "wt")
    else:
        fw = open(outfile, "w+")
    
    alignments = load_alignments(infile1, chrom, bin_start, bin_end)
    bundles = BundleBuilder(alignments, keep=True)
    for n, bundle in enumerate(bundles):
        # if n % 100 == 0:
        #     pass
            # print(n)
        start = bundle.start_min
        end = bundle.end_max
        width = end - start
        data = defaultdict(list)
        for obj in bundle.data:
            rg = obj.segment.get_tag("RG")
            data[rg].append(obj)
        array1 = [] # Read groups
        array2 = [] # alignments
        for k, v in data.items():
            array1.append(k)
            array2.append(v)
        array3 = [] # coverages
        for objs in array2:
            covs = np.zeros(width, dtype=np.int)
            for obj in objs:
                for idx in np.arange(obj.start - start, obj.end - start):
                    covs[idx] += 1
            array3.append(covs)
        array4 = [] # events
        for objs in array2:
            events = defaultdict(list)
            for obj in objs:
                for e in obj.me:
                    if e[1] == "-":
                        continue
                    for bi in range(len(e[1])):
                        events[e[0] + bi].append(e[2])
            array4.append(events)
        for pos in range(start, end):
            if pos < bin_start:
                continue
            elif pos >= bin_end:
                break
            # print(pos, start0, end0)
            idx = pos - start
            covs_all = [item[idx] for item in array3]
            covs = []
            cell_count = 0
            read_count = 0
            ref = sequence[pos - bin_start]
            array5 = [] # base counter
            scores = []
            for i, cov in enumerate(covs_all):
                if cov == 0:
                    continue
                covs.append(cov)
                cell_count += 1
                read_count += cov
                ct = Counter(array4[i][pos])
                tot = sum(ct.values())
                ref_count = cov - tot
                if ref_count > 0:
                    ct[ref] = ref_count
                array5.append(ct)
            
            posible = []
            for ct in array5:
                for k in ct.keys():
                    posible.append(k)
            posible = list(set(posible))
            
            counter = dict()
            for base in posible:
                count = 0
                for ct in array5:
                    count += ct[base]
                counter[base] = count
            allele = None
            allele_count = None
            count_max = max(counter.values())
            ambiguous = False
            for k, v in counter.items():
                if allele is None:
                    allele = k
                    allele_count = v
                else:
                    if v > allele_count:
                        allele = k
                        allele_count = v
                    elif v == allele_count and v == count_max:
                        ambiguous = True
                        break
            score = 0
            for cov, ct in zip(covs, array5):
                c = ct[allele]
                r = c / cov
                if c >= 2 and r >= 0.8:
                    score += 1

            conf = cell_count >= 2 and score == cell_count
            
            s1 = "."
            s2 = "."
            # if allele_count != read_count:
            if True:
                items1 = []
                for ct in array5:
                    items2 = []
                    for k, v in ct.items():
                        if v == 0:
                            continue
                        items2.append("%s:%d" % (k, v))
                    items1.append(",".join(items2))
                s1 = ";".join(items1)

                items1 = []
                for k, v in counter.items():
                    items1.append("%s:%d" % (k, v))
                s2 = ",".join(items1)
            conf = "T" if conf else "F"
            line = "\t".join(map(str, [pos, ref, allele, allele_count, cell_count, read_count, score, conf, s2, s1]))
            fw.write(line + "\n")
        # break
    fw.close()  
        
    # print(task)

def main():
    # input.bam genome.fasta chr1 8 outdir
    infile1, infile2, chrom, threads, outdir = sys.argv[1:]
    threads = int(threads)
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    with pysam.AlignmentFile(infile1) as f:
        length = f.get_reference_length(chrom)
        
    tasks = []
    width = 1000000 # bin width
    for i, start in enumerate(range(0, length, width)):
        if start >= length:
            break
        end = min(start + width, length)
        outfile = os.path.join(outdir, "%04d.matrix.gz" % i)
        task = {
            "bam": infile1,
            "fasta": infile2,
            "chrom": chrom,
            "start": start,
            "end": end,
            "outfile": outfile}
        tasks.append(task)
    
    with open(os.path.join(outdir, "meta.tsv"), "w+") as fw:
        for task in tasks:
            line = "\t".join(map(str, [
                task["chrom"], task["start"], task["end"], task["outfile"]]))
            fw.write(line + "\n")
                
    if threads == 1:
        for task in tasks:
            execute_bin(task)
    else:
        pool = multiprocessing.Pool(threads)
        rets = [pool.apply_async(execute_bin, (t,)) for t in tasks]
        pool.close()
        pool.join()
        for ret in rets:
            assert ret.successful() == True
    
    
if __name__ == '__main__':
    main()
    

      