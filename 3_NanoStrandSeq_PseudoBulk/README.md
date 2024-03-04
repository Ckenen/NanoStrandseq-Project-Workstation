# Analysis of NanoStrand-seq as pseudobulk

## Preparation

Download PB-CCS and ONT-UL BAMs for HG001 from GIAB:

    mkdir -p data
    cd data
    wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/PacBio_SequelII_CCS_11kb/HG001_GRCh38/HG001_GRCh38.haplotag.RTG.trio.bam
    wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/PacBio_SequelII_CCS_11kb/HG001_GRCh38/HG001_GRCh38.haplotag.RTG.trio.bam.bai
    wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Ultralong_OxfordNanopore/NA12878-minion-ul_GRCh38.bam
    wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Ultralong_OxfordNanopore/NA12878-minion-ul_GRCh38.bam.bai
    cd ..

The benchmark variant calls is from https://github.com/Ckenen/GRCh38-HG001-Variant-Calls

Download genome stratification files for GRCh38 from GIAB:

    cd data
    wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/v3.1-genome-stratifications-GRCh38.tar.gz
    tar -zxvf v3.1-genome-stratifications-GRCh38.tar.gz
    cd ..

Deposited NanoStrand-seq BAM files in data/nss_bams directory.

    mkdir -p data/nss_bams
    cd data/nss_bams
    # 
    cd ..

## 

