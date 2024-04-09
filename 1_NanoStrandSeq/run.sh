#!/bin/sh

for smk in 1_SnakeQC.smk 2_SnakeDemux.smk 3_SnakeMapping.smk 4_SnakeCount.smk 5_SnakeStat.smk 6_SnakePreSeq.smk; do
    echo "Running: $smk, configfile: $f"

    # Test
    snakemake -s ${smk} -np

    # Running at cluster
    # ./${smk}

    # Running at local
    # snakemake -s ${smk} --use-conda -j

done
