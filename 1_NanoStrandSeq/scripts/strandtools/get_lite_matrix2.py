#!/usr/bin/env python
import sys
from collections import defaultdict


def get_base_counter(s):
    base_counter = defaultdict(int)
    if s != ".":
        for cell_s in s.split(";"):
            d_cell = dict()
            for base_count in cell_s.split(","):
                base, count = base_count.split(":")
                count = int(count)
                d_cell[base] = count
                base_counter[base] += count
    return base_counter


def get_potential_base(base_counter):
    base = None
    if len(base_counter) > 0:
        items = list(sorted(base_counter.items(), key=lambda item: item[1]))
        base = items[-1][0]
    return base


def main():
    for line in sys.stdin:
        row = line.strip("\n").split("\t")

        ref_base = row[1]
        
        if row[7] == "." or row[8] == ".":
            continue

        base_counter_1 = get_base_counter(row[7])
        base_counter_2 = get_base_counter(row[8])
        base_depth_1 = sum(base_counter_1.values())
        base_depth_2 = sum(base_counter_2.values())

        assert base_depth_1 > 0 and base_depth_2 > 0

        potential_base_1 = get_potential_base(base_counter_1)
        potential_base_2 = get_potential_base(base_counter_2)

        if potential_base_1 == ref_base \
            and potential_base_2 == ref_base \
            and base_counter_1[potential_base_1] >= base_depth_1 * 0.8 \
            and base_counter_2[potential_base_2] >= base_depth_2 * 0.8:
            continue
            
        sys.stdout.write(line)

if __name__ == "__main__":
    main()