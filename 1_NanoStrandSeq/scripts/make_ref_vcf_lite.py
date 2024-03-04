#!/usr/bin/env python
import sys
import re


def main():

    for line in sys.stdin:
        line = line.strip("\n")
        if line.startswith("#"):
            if line.startswith("##contig"):
                if re.search("ID=chr([0-9]+|[XYM]),", line):
                    print(line)
            elif line.startswith("##FILTER"):
                continue
            elif line.startswith("##INFO"):
                continue
            else:
                print(line)
        else:
            row = line.split("\t")
            if re.search("^chr([0-9]+|[XYM])$", row[0]) is None:
                continue
            row[6] = "."
            row[7] = "."
            line = "\t".join(row)
            print(line)


if __name__ == '__main__':
    main()
