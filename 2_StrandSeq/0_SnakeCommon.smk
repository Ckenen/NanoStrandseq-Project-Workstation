#!/usr/bin/env
import pandas as pd
configfile: "config.yaml" 
THREADS = config["THREADS"]
RUNS = config["RUNS"]
DAT = pd.read_excel(config["TABLE"])
DAT.index = DAT["Cell"]
DAT = DAT[DAT["Run"].isin(RUNS)]
DAT["RunCell"] = ["%s/%s" % (run, cell) for run, cell in DAT[["Run", "Cell"]].values]
RUN_CELLS = list(DAT["RunCell"])
print("runs: %d" % len(RUNS))
print("Cells: %d" % len(RUN_CELLS))