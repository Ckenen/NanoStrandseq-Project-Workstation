#!/usr/bin/env python
import sys
from collections import defaultdict
import pysam


def load_segments(path):
    items = []
    with pysam.AlignmentFile(path) as f:
        for segment in f:
            item = [segment.query_name, segment.mapping_quality, segment.reference_name]
            items.append(item)
    items = list(sorted(items, key=lambda item: item[0]))
    return items


def main():
    infile, outfile = sys.argv[1:]
    
    counter = defaultdict(int)

    items = None
    for item in load_segments(infile):
        if items is None:
            items = [item]
        elif item[0] == items[0][0]:
            items.append(item)
        elif item[0] > items[0][0]:
            if len(items) > 1:
                items = list(sorted(items, key=lambda item: item[1]))
            counter[items[-1][2].split("_")[0]] += 1
            items = [item]
        else:
            assert False
    if items is not None:
        if len(items) > 1:
            items = list(sorted(items, key=lambda item: item[1]))
        counter[items[-1][2].split("_")[0]] += 1
            
    with open(outfile, "w+") as fw:
        s = sum(counter.values())
        for k, v in sorted(counter.items()):
            if s == 0:
                r = 0
            else:
                r = v / s
            fw.write("%s\t%d\t%f\n" % (k, v, r))
    
    
if __name__ == '__main__':
    main()
    