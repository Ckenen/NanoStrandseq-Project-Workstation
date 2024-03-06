#!/bin/sh

for f in config_HG001.yaml config_C57DBA.yaml; do
    #./1_SnakeQC.smk --configfile $f
    #./2_SnakeDemux.smk --configfile $f
    #./3_SnakeMapping.smk --configfile $f
    #./4_SnakeCount.smk --configfile $f
    #./5_SnakeStat.smk --configfile $f
    ./6_SnakePreSeq.smk --configfile $f
done


## Example of local running:
# snakemake -s 1_SnakeQC.smk --configfile config_HG001.yaml -j
