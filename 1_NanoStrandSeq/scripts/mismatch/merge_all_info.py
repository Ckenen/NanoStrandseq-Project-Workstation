#!/usr/bin/env python
import sys
import gzip


def load_mismatch_sites(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            yield line.strip("\n")


def load_base_coverages(paths):
    handles = [gzip.open(p, "rt") for p in paths]
    for items in zip(*handles):
        yield ":".join([line.strip("\n") for line in items])
    for h in handles:
        h.close()


def load_mismatch_events(paths):
    handles = [gzip.open(p, "rt") for p in paths]
    for items in zip(*handles):
        yield ":".join([line.strip("\n") for line in items])
    for h in handles:
        h.close()


def main():
    in_tsv = sys.argv[1]
    out_tsv = sys.argv[-1]
    paths = sys.argv[2:-1]
    n = int(len(paths) / 2)

    cov_tsv_list = paths[:n]
    event_tsv_list = paths[n:]
    loader1 = load_mismatch_sites(in_tsv)
    loader2 = load_base_coverages(cov_tsv_list)
    loader3 = load_mismatch_events(event_tsv_list)
    with gzip.open(out_tsv, "wt") as fw:
        for values in zip(loader1, loader2, loader3):
            line = "\t".join(values)
            fw.write(line)
            fw.write("\n")


if __name__ == "__main__":
    main()
