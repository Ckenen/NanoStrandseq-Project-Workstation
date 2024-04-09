#!/bin/sh

for smk in 1_SnakePrepare.smk 2_SnakeMapping.smk 3_SnakeCount.smk 4_SnakeStat.smk 5_SnakePreSeq.smk; do
    echo $smk

    # Test
    snakemake -s $smk -np

    # Local
    # snakemake -s $smk -j --use-conda

    # Cluster
    # ./${smk}

done
