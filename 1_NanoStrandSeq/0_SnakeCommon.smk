#!/usr/bin/env runsnakemake
import pandas as pd
configfile: "config.yaml" 
threads = config["threads"]
runs = config["runs"]

# Uniq reads
path = "results/mapping/final_reads.%s.tsv" % config["name"]
cell2reads = dict()
if os.path.exists(path):
    tmp = pd.read_csv(path, sep="\t")
    cell2reads = {cell: reads for cell, reads in tmp[["Cell", "UniqReads"]].values}

# Table
dat = pd.read_excel(config["table"])
dat.index = dat["Cell"]
dat = dat[[run in runs for run in dat["Run"]]]
dat["RunCell"] = ["%s/%s" % (run, cell) for run, cell in dat[["Run", "Cell"]].values]
dat["Reads"] = [cell2reads.get(cell) for cell in dat["Cell"]]
run_cells = list(dat["RunCell"])
run_cells_filtered = list(dat[~(dat["Reads"] < config["min_reads"])]["RunCell"])
print("Runs: %d" % len(runs))
print("Cells: %d" % len(run_cells))
print("Cells (fitlered): %d" % len(run_cells_filtered))