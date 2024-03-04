#!/usr/bin/env python
import sys


def get_reads(data):
    return sum([item[1] for item in data])


def get_total_length(data):
    return sum([length * count for length, count in data])


def get_mean(data):
    bases = 0
    reads = 0
    for length, count in data:
        bases += length * count
        reads += count
    mean = bases / reads
    return mean

    
def get_median(data):
    reads = get_reads(data)
    median = None
    if reads % 2 == 1:
        i = int(reads / 2)
        j = 0
        for length, count in data:
            j += count
            if j > i:
                median = length
                break
    else:
        i2 = int(reads / 2)
        i1 = i2 - 1
        j = 0
        v1 = None
        for length, count in data:
            j += count
            if j > i1 and v1 is None:
                v1 = length
            if j > i2:
                median = (length + v1) / 2
                break
    return median
    
    
def get_n50(data):
    n50 = None
    i = get_total_length(data) / 2
    j = 0
    for length, count in data[::-1]:
        j += length * count
        if j >= i:
            n50 = length
            break
    return n50

    
def main():
    infile = sys.argv[1]
    
    run = infile.split("/")[-1][:-len(".txt")]

    data = []
    with open(infile) as f:
        for line in f:
            length, count = line.strip("\n").split("\t")
            data.append([int(length), int(count)])
    data.sort()

    total_length = get_total_length(data)
    reads = get_reads(data)
    mean = get_mean(data)
    median = get_median(data)
    n50 = get_n50(data)
    
    print("Run", "TotalBase", "TotalRead", "MaxLength",
          "MeanLength", "MedianLength", "N50", sep="\t")
    print(run, total_length, reads, data[-1][0], mean, median, n50, sep="\t")


if __name__ == '__main__':
    main()
