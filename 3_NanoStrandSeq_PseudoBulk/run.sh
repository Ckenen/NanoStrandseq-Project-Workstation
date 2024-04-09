#!/usr/bin/env bash

for smk in 1_SnakeMakeBam.smk 2_SnakeStat.smk 3_SnakeSNV.smk 4_SnakeSV.smk 5_SnakeStratification.smk; do
    snakemake -s ${smk} -j --use-conda
done