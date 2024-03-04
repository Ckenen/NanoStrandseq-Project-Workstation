# Test Whatshap

    mkdir -p test_whatshap
    
    /date/chenzonggui/baixiuzhen/GRCh38-HG001-Variant-Calls/results/benchmark_autosomal_v4.2.1_chrx_v3.3.2.vcf.gz

    /date/chenzonggui/baixiuzhen/GIAB/HG001_GRCh38_NISTv4.2.1/PacBio_SequelII_CCS_11kb/HG001_GRCh38.haplotag.RTG.trio.bam

    /datb1/chenzonggui/NanoStrandSeq_Sup/assembly/HG001_Cell300/round2/snvs.vcf.gz

    whatshap unphase /date/chenzonggui/baixiuzhen/GRCh38-HG001-Variant-Calls/results/benchmark_autosomal_v4.2.1_chrx_v3.3.2.vcf.gz | bgzip -c > test_whatshap/unphased.vcf.gz
    tabix -p vcf test_whatshap/unphased.vcf.gz

    whatshap phase -r /home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa test_whatshap/unphased.vcf.gz /date/chenzonggui/baixiuzhen/GIAB/HG001_GRCh38_NISTv4.2.1/PacBio_SequelII_CCS_11kb/HG001_GRCh38.haplotag.RTG.trio.bam /datb1/chenzonggui/NanoStrandSeq_Sup/assembly/HG001_Cell300/round2/snvs.vcf.gz | bgzip -c > test_whatshap/phased.vcf.gz

    whatshap phase -r /home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.canonical.genome.fa --chromosome chr1 --ignore-read-groups test_whatshap/unphased.vcf.gz /datb1/chenzonggui/NanoStrandSeq_Sup/assembly/HG001_Cell300/round2/snvs.vcf.gz /date/chenzonggui/baixiuzhen/GIAB/HG001_GRCh38_NISTv4.2.1/PacBio_SequelII_CCS_11kb/HG001_GRCh38.haplotag.RTG.trio.bam | bgzip -c > test_whatshap/phased.chr1.v2.vcf.gz


set +u; source activate sniffles2
        sniffles -t {threads} --phase --output-rnames --tandem-repeats {input.bed} \
            --reference {input.fsa} -i {input.bam} -v {output.vcf1}
        cat {output.vcf1} | grep '#' > {output.vcf2}
        cat {output.vcf1} | grep -v '#' | grep -v 'INV' | grep -v 'DUP' | grep -v 'BND' >> {output.vcf2}
        bgzip -c {output.vcf2} > {output.vcf3}
        tabix -p vcf {output.vcf3} ) &> {log}

sniffles -t 24 --phase --output-rnames --tandem-repeats /home/chenzonggui/species/homo_sapiens/TandemRepeats/GRCh38.TandemRepeats.bed --reference /home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.primary_assembly.genome.fa -i results/HG001_Cell_350/prepare/all_cells.all_chroms.raw.bam -v HG001_Cell_350.sniffles2.vcf

cat HG001_Cell_350.sniffles2.vcf | grep '#' > HG001_Cell_350.sniffles2.filtered.vcf
cat HG001_Cell_350.sniffles2.vcf | grep -v '#' | grep -v 'INV' | grep -v 'DUP' | grep -v 'BND' >> HG001_Cell_350.sniffles2.filtered.vcf
bgzip HG001_Cell_350.sniffles2.filtered.vcf
tabix -p vcf HG001_Cell_350.sniffles2.filtered.vcf.gz

sniffles -t 24 --phase --output-rnames --genotype-vcf HG001_Cell_350.sniffles2.filtered.vcf.gz --tandem-repeats /home/chenzonggui/species/homo_sapiens/TandemRepeats/GRCh38.TandemRepeats.bed --reference /home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.primary_assembly.genome.fa -i results/HG001_Cell_350/round2/merged.hp1.bam -v HG001_Cell_350.sniffles2.hp1.vcf

cat HG001_Cell_350.sniffles2.hp1.vcf | grep '#' > HG001_Cell_350.sniffles2.hp1.filtered.vcf
cat HG001_Cell_350.sniffles2.hp1.vcf | grep -v '#' | grep -v 'INV' | grep -v 'DUP' | grep -v 'BND' >> HG001_Cell_350.sniffles2.hp1.filtered.vcf
bgzip HG001_Cell_350.sniffles2.hp1.filtered.vcf
tabix -p vcf HG001_Cell_350.sniffles2.hp1.filtered.vcf.gz

sniffles -t 24 --phase --output-rnames --genotype-vcf HG001_Cell_350.sniffles2.filtered.vcf.gz --tandem-repeats /home/chenzonggui/species/homo_sapiens/TandemRepeats/GRCh38.TandemRepeats.bed --reference /home/chenzonggui/species/homo_sapiens/GRCh38.p13/GRCh38.primary_assembly.genome.fa -i results/HG001_Cell_350/round2/merged.hp2.bam -v HG001_Cell_350.sniffles2.hp2.vcf

cat HG001_Cell_350.sniffles2.hp2.vcf | grep '#' > HG001_Cell_350.sniffles2.hp2.filtered.vcf
cat HG001_Cell_350.sniffles2.hp2.vcf | grep -v '#' | grep -v 'INV' | grep -v 'DUP' | grep -v 'BND' >> HG001_Cell_350.sniffles2.hp2.filtered.vcf
bgzip HG001_Cell_350.sniffles2.hp2.filtered.vcf
tabix -p vcf HG001_Cell_350.sniffles2.hp2.filtered.vcf.gz
