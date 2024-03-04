#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

f_tsv, f_pdf = sys.argv[1:]

dat = pd.read_csv(f_tsv, sep="\t")
nbin = dat["Bin"].max() + 1
vs = np.concatenate([dat["Crick"], dat["Watson"]])
mean = np.mean(vs)
std = np.std(vs)
vmax = max(dat["Crick"].max(), dat["Watson"].max())
xlim = vmax * 1.2
xlim = mean + 3 * std
xlim = max(10, xlim)

clusters = list(sorted(set(dat["Name"]), key=lambda item: int(item.split("_")[1])))

ncol = 10
nrow = int(len(clusters) / ncol)
if len(clusters) % ncol > 0:
    nrow += 1

fig, axs = plt.subplots(nrow, ncol, figsize=(11, nrow * 2.5), sharex=True, sharey=True)

for i, cluster in enumerate(clusters):
    i1 = int(i / ncol)
    i2 = i % ncol
    ax = axs[i1][i2]
    plt.sca(ax)
    
    plt.xlabel(cluster)
    
    tmp = dat[dat["Name"] == cluster]
    xs1 = tmp["Crick"]
    xs2 = tmp["Watson"]
    ys = np.arange(len(tmp))
    
    plt.barh(ys, xs1, color="C1", height=1)
    plt.barh(ys, -xs2, color="C0", height=1)
    plt.plot([0, 0], [0, len(ys)], lw=2, color="grey")
    plt.xlim(-xlim, xlim)
    plt.ylim(-0.5, nbin - 0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # break
    
plt.tight_layout()
plt.savefig(f_pdf, dpi=300)