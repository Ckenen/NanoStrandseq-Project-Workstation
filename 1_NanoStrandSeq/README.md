# Snakemake workflow for basic analysis of NanoStrand-seq

Commands:

    snakemake -s 1_SnakeQC.smk -np
    snakemake -s 2_SnakeDemux.smk -np
    snakemake -s 3_SnakeMapping.smk -np
    snakemake -s 4_SnakeBarPlot.smk -np

samtools sort -@ 4 -n -o sorted_by_name.bam results/mapping/minimap2/20220708_GM12878/20220708_GM12878.sc001.bam

sstools FilterBam -n '^chr([0-9]+|[XY])$' -q 30 results/mapping/minimap2/20220708_GM12878/20220708_GM12878.sc001.bam 20220708_GM12878.sc001_filtered.bam 
samtools index -@ 4 20220708_GM12878.sc001_filtered.bam