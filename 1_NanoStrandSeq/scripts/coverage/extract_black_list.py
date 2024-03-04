#!/usr/bin/env python
import sys
import numpy as np
import pysam


MAX_COV = 10    # depth超过10的区间为感兴趣的区间
MAX_GAP = 100   # 分批处理时，批次之间的坐标的最小间隔
MAX_MERGE_GAP = 50 # 兴趣区间小于该值将会被合并起来

def process(array):
    if len(array) < MAX_COV:
        return None
    start = array[0].reference_start
    end = max([segment.reference_end for segment in array])
    values = np.zeros(end - start, dtype=int)
    for segment in array:
        for i in range(segment.reference_start - start, segment.reference_end - start):
            values[i] += 1
    if max(values) < MAX_COV:
        return None
    intervals = []
    begin = None
    for i, v in enumerate(values):
        if v >= MAX_COV:
            if begin is None:
                begin = i
        else:
            if begin:
                intervals.append([begin, i])
                begin = None
    if begin:
        intervals.append([begin, len(values)])
    intervals = [[s + start, e + start] for s, e in intervals]
    i = 1
    while i < len(intervals):
        if intervals[i][0] - intervals[i - 1][1] < MAX_MERGE_GAP:
            intervals[i - 1][1] = intervals[i][1]
            intervals.pop(i)
        else:
            i += 1
    for s, e in intervals:
        print(array[0].reference_name, s, e, "Blacklist", ".", "+", sep="\t")


def main():
    infile = sys.argv[1]
    array = None
    max_end = 0
    with pysam.AlignmentFile(infile) as f:
        for segment in f:
            if array:
                if segment.reference_name == array[0].reference_name:
                    if segment.reference_start >= max_end + MAX_GAP:
                        process(array)
                        array = [segment]
                        max_end = segment.reference_end
                    else:
                        array.append(segment)
                        max_end = max(max_end, segment.reference_end)
                else:
                    process(array)
                    array = [segment]
                    max_end = segment.reference_end
            else:
                array = [segment]
                max_end = segment.reference_end
    if array:
        process(array)


if __name__ == '__main__':
    main()
