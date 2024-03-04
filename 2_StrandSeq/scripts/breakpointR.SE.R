#!/usr/bin/env Rscript
library(breakpointR)
args <- commandArgs(trailingOnly = TRUE)
datafolder <- args[1]
threads <-  args[2]
outputfolder <- args[3]

threads <- as.numeric(threads)

## Run breakpointR
breakpointr(inputfolder = datafolder, 
    outputfolder = outputfolder,
    pairedEndReads = FALSE, 
    numCPU=threads,
    reuse.existing.files = FALSE, 
    windowsize = 20,
    binMethod = 'reads', 
    pair2frgm = FALSE, 
    min.mapq = 1,
    filtAlt = TRUE)
