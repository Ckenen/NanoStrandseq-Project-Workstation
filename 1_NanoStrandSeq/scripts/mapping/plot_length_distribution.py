#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")


def main():
    infile, outfile = sys.argv[1:]

    with open(infile) as f:
        values = [int(line.strip("\n")) for line in f]
    mean = np.mean(values)
    median = np.median(values)
    std = np.std(values)

    WIDTH = 20
    MAX_VALUE = 5000

    xs = np.arange(0, MAX_VALUE, WIDTH) + int(WIDTH / 2)
    ys = np.zeros(len(xs))
    for v in values:
        i = int(v / WIDTH)
        if 0 <= i < len(ys):
            ys[i] += 1
    max_y = max(ys)

    plt.figure(figsize=(5, 3))
    plt.bar(xs, ys, width=WIDTH)
    plt.text(MAX_VALUE * 0.65, max_y * 0.9, "N=%d" % len(values))
    plt.text(MAX_VALUE * 0.65, max_y * 0.8, "Mean=%.2f" % mean)
    plt.text(MAX_VALUE * 0.65, max_y * 0.7, "Median=%d" % median)
    plt.text(MAX_VALUE * 0.65, max_y * 0.6, "Std=%.2f" % std)
    plt.xlabel("Read length (nt)")
    plt.ylabel("Read number")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
