#!/usr/bin/env runsnakemake
import pandas as pd
configfile: "config.yaml" 
threads = config["threads"]
runs = config["runs"]
dat = pd.read_excel(config["table"])
dat.index = dat["Cell"]
dat = dat[[run in runs for run in dat["Run"]]]
dat["RunCell"] = ["%s/%s" % (run, cell) for run, cell in dat[["Run", "Cell"]].values]
run_cells = list(dat["RunCell"])
print("Runs: %d" % len(runs))
print("Cells: %d" % len(run_cells))
