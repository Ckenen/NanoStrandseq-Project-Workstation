#!/bin/sh

bam=$1
fasta=$2
chrom=$3
threads=$4
vcf=$5

dirname=`dirname $vcf`
basename=`basename $vcf .vcf.gz`
vcf_tmp="${dirname}/${basename}.vcf"

ont_model="/home/chenzonggui/software/princess/bin/modules/ont/model"
model="/home/chenzonggui/software/princess/bin/modules/ont/model"

set +u
source activate clair

clair.py callVarBam --delay 0 --chkpnt_fn ${model} --ref_fn ${fasta} --bam_fn ${bam} --ctgName ${chrom} --sampleName ${basename} --threads ${threads} --call_fn ${vcf_tmp}

conda deactivate

bgzip -c ${vcf_tmp} > ${vcf}
tabix -p vcf ${vcf}
