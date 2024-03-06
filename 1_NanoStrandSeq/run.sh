#!/bin/sh

for f in config_HG001.yaml config_B6D2F1.yaml; do
    for smk in 1_SnakeQC.smk 2_SnakeDemux.smk 3_SnakeMapping.smk 4_SnakeCount.smk 5_SnakeStat.smk 6_SnakePreSeq.smk; do
        echo "Running: $smk, configfile: $f"

        # Test
        snakemake -s ${smk} --configfile $f -np

        # Running at cluster
        # ./${smk} --configfile $f

        # Running at local
        # snakemake -s ${smk} --configfile $f -j

    done
done