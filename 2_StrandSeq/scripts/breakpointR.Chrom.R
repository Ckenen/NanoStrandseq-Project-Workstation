#!/usr/bin/env Rscript

library(breakpointR)

args <- commandArgs(trailingOnly = TRUE)

datafolder <- args[1]
chromosome <- args[2]
threads <-  args[3]
outputfolder <- args[4]

threads <- as.numeric(threads)

## Run breakpointR
breakpointr(inputfolder = datafolder, 
    outputfolder = outputfolder,
    chromosomes = chromosome, 
    pairedEndReads = FALSE, 
    numCPU=threads,
    reuse.existing.files = FALSE, 
    windowsize = 20,
    binMethod = 'reads', 
    pair2frgm = FALSE, 
    min.mapq = 1,
    filtAlt = TRUE)
