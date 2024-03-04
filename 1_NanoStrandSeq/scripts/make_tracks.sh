#!/bin/sh

# prefix.bed.gz
# prefix.bed.gz.tbi
# prefix.bw
# prefix.+.bw
# prefix.-.bw

bam=$1
prefix=$2


tmpdir="${prefix}.sort_tmp"
mkdir -p ${tmpdir}

bed="${prefix}.bed"
bed_gz="${prefix}.bed.gz"
bedtools bamtobed -bed12 -i ${bam} \
    | awk -v FS='\t' -v OFS='\t' '{if($6=="+"){c="107,137,138"}else{c="248,173,97"}{print $1,$2,$3,$4,$5,$6,$7,$8,c,$10,$11,$12}}' \
    | sort -k1,1 -k2,2n -k3,3n -T ${tmpdir} > ${bed}
bgzip -c ${bed} > ${bed_gz}
tabix -p bed ${bed_gz}

txt="${prefix}.txt"
samtools view -H ${bam} | grep '@SQ' | awk -v OFS='\t' '{print $2,$3}' | sed 's/SN://g' | sed 's/LN://g' > ${txt}
        
bg="${prefix}.bg"
bg1="${prefix}.+.bg"
bg2="${prefix}.-.bg"

bedtools genomecov -bga -i ${bed} -g ${txt} | sort -k1,1 -k2,2n -k3,3n -T ${tmpdir} > ${bg}
bedtools genomecov -bga -strand + -i ${bed} -g ${txt} | sort -k1,1 -k2,2n -k3,3n -T ${tmpdir} > ${bg1}
bedtools genomecov -bga -strand - -i ${bed} -g ${txt} | sort -k1,1 -k2,2n -k3,3n -T ${tmpdir} > ${bg2}

bw="${prefix}.bw"
bw1="${prefix}.+.bw"
bw2="${prefix}.-.bw"

bedGraphToBigWig ${bg} ${txt} ${bw}
bedGraphToBigWig ${bg1} ${txt} ${bw1}
bedGraphToBigWig ${bg2} ${txt} ${bw2}

rm ${bed} ${txt} ${bg} ${bg1} ${bg2}
rm -r ${tmpdir}