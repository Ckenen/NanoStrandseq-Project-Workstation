#!/usr/bin/env python
import sys
import multiprocessing
import pysam


def process_chromosome(bamfile, chrom, width, step):
    with pysam.AlignmentFile(bamfile) as f:
        length = f.get_reference_length(chrom)
        num = int((length - width) / step) + 1
        if (length - width) % step > 0:
            num += 1
        bins = []
        for i in range(num):
            start = step * i
            end = min(start + width, length)
            bins.append([chrom, start, end, "Bin.%d" % i, 0, "+"])
        for n, segment in enumerate(f.fetch(chrom)):
            # if n >= 10000:
            #     break
            start = segment.reference_start
            end = segment.reference_end
            i1 = int(start / step)
            i2 = int((start - width) / step) + 1
            i3 = int((end - 1) / step)
            i4 = int((end - 1 - width) / step) + 1
            i5, i6 = min(i1, i2), max(i3, i4) + 1
            assert i5 < i6
            for i in range(i5, i6):
                bins[i][4] += 1
    return bins


def main():
    bamfile, threads, outfile = sys.argv[1:]
    threads = int(threads)
    width = 10000 # 10k
    step = 5000

    results = []
    pool = multiprocessing.Pool(threads)
    with pysam.AlignmentFile(bamfile) as f:
        for chrom in f.references:
            args = (bamfile, chrom, width, step)
            r = pool.apply_async(process_chromosome, args)
            results.append(r)
    pool.close()
    pool.join()

    results = [r.get() for r in results]

    with open(outfile, "w+") as fw:
        for rows in results:
            for row in rows:
                line = "\t".join(map(str, row))
                fw.write(line + "\n")


if __name__ == "__main__":
    main()