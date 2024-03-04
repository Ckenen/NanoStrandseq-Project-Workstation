#!/usr/bin/env Rscript
library(breakpointR)
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]
load(infile)
write.table(breakpoints$lib.metrics, outfile, sep = "\t")