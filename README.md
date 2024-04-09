# NanoStrandseq-Project-Workstation

This repository includes the source code to plot figures in the paper and the analysis workflows of NanoStrand-seq and Strand-seq (as described below).

| Directory | Description |
| :-- | :-- |
| 0_PlotFigures | Source codes to plot figures in paper. |
| 1_NanoStrandSeq | Workflow of basic analysis of NanoStrand-seq datasets. |
| 2_StrandSeq | Workflow of basic analysis of Strand-seq datasets. |
| 3_NanoStrandSeq_PseudoBulk | Workflow of pseudo-bulk analysis of NanoStrand-seq datasets. |
| 4_NanoStrandSeq_Phasing | Workflow of phasing SNPs and SVs by NanoStrand-seq. |
| 5_NanoStrandSeq_Assembly | Workflow of de novo assembly by HiFi and NanoStrand-seq. |


## Dependencies

Benchmark small variant calls (SNPs and Indels) and regions for HG001:

https://github.com/Ckenen/GRCh38_HG001_SNP_Indel

Benchmark structure variant calls (insertions and deletions) and regions for HG001:

https://github.com/Ckenen/GRCh38_HG001_SV

Benchmark small variant calls (SNPs and Indels) and regions for B6D2F1 mouse:

https://github.com/Ckenen/GRCm38_B6D2F1_SNP_Indel

Benchmark structure variant calls (insertions and deletions) and regions for B6D2F1 mouse:

https://github.com/Ckenen/GRCm38_B6D2F1_SV

Find barcodes in lone-reads:

https://github.com/Ckenen/fbilr

Split and trim single-cell reads:

https://github.com/Ckenen/nss-demultiplexing

Toolkits for Strand-seq and NanoStrand-seq:

https://github.com/Ckenen/sstools

Others dependencies:

https://github.com/Ckenen/pyBioInfo

