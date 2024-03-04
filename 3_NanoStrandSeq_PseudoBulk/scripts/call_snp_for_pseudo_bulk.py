#!/usr/bin/env python
import optparse
import multiprocessing as mp
from collections import defaultdict
import pysam
from pyBioInfo.Range import GRange
from pyBioInfo.IO.File import Alignment
from pyBioInfo.Utils import ShiftLoader


def make_vcf_grange(record):
    obj = GRange(chrom=record.chrom, start=record.start, end=record.start + 1)
    obj.record = record
    return obj

     
def call_snps(f_bam, f_vcf, chrom, start, end):
    new_snps = []
    with pysam.VariantFile(f_vcf) as vcf, pysam.AlignmentFile(f_bam) as bam:
        snps = [make_vcf_grange(s) for s in vcf.fetch(chrom, start, end)]
        snps = list(filter(lambda snp: len(snp.record.ref) == 1, snps))
        alignments = [Alignment(s) for s in bam.fetch(chrom, start, end)]
        loader = ShiftLoader(alignments)
        for snp in snps:
            base_counter = defaultdict(list)
            for alignment in loader.fetch(obj=snp):
                rg = alignment.segment.get_tag("RG")
                base = alignment.get_query_base(snp.start)
                base_counter[base].append(rg)
            ref_base = snp.record.ref
            base_reads = {base: len(rgs) for base, rgs in base_counter.items()}
            base_cells = {base: len(set(rgs)) for base, rgs in base_counter.items()}
            items = list(sorted(base_reads.items(), key=lambda item: item[1], reverse=True))
            total_reads = sum(base_reads.values())
            
            new_snp = dict()
            new_snp["chrom"] = chrom
            new_snp["start"] = snp.start
            new_snp["ref_base"] = ref_base
            new_snp["id"] = snp.record.id
            new_snp["base_counter"] = base_counter
            new_snp["base_reads"] = base_reads
            new_snp["base_cells"] = base_cells
            
            if len(items) == 0:
                continue   
            elif len(items) == 1:
                if items[0][0] == ref_base:
                    continue
                if items[0][0] == "-":
                    continue
                base, reads = items[0]
                cells = base_cells[base]
                if reads >= 4 and cells >= 2:
                    new_snp["base1"] = base
                    new_snp["base2"] = base
                    new_snp["genotype"] = "HOM"
                else:
                    continue
            else:
                base1, reads1 = items[0]
                base2, reads2 = items[1]
                cells1 = base_cells[base1]
                cells2 = base_cells[base2]
                ratio1 = reads1 / total_reads
                ratio2 = reads2 / total_reads
                if base1 == "-":
                    if base2 == ref_base:
                        continue
                    else:
                        if reads2 >= 4 and cells2 >= 2 and ratio2 >= 0.66:
                            new_snp["base1"] = base2
                            new_snp["base2"] = base2
                        else:
                            continue
                elif base2 == "-":
                    if base1 == ref_base:
                        continue
                    else:
                        if reads1 >= 4 and cells1 >= 2 and ratio1 >= 0.66:
                            new_snp["base1"] = base1
                            new_snp["base2"] = base1
                        else:
                            continue
                else:
                    if reads1 >= 4 and cells1 >= 2 and ratio1 >= 0.8:
                        if base1 == ref_base:
                            continue
                        else:
                            new_snp["base1"] = base1
                            new_snp["base2"] = base1
                    elif reads1 >= 4 and cells1 >= 2 and reads2 >= 4 and cells2 >= 2 and ratio1 >= 0.4 and ratio2 >= 0.4:
                        new_snp["base1"] = base1
                        new_snp["base2"] = base2
                    else:
                        continue
            new_snps.append(new_snp)
    return new_snps


def main():
    parser = optparse.OptionParser(usage="%prog [options] input.bam")
    
    parser.add_option("-o", "--output", dest="outvcf")
    parser.add_option("-s", "--snps", dest="snpvcf")
    parser.add_option("-t", "--threads", dest="threads", default=1, type="int")
    parser.add_option("-c", "--chrom", dest="chrom")
    
    options, args = parser.parse_args()
    
    f_vcf = options.snpvcf
    f_bam = args[0]
    f_out = options.outvcf
    chrom = options.chrom
    threads = options.threads
    width = 1000000
    sample = "RESULTS"
    
    assert f_vcf
    
    task_list = []
    header = None
    with pysam.AlignmentFile(f_bam) as bam:
        header = bam.header.as_dict()
        if chrom is None:
            for c in sorted(bam.references):
                length = bam.get_reference_length(c)
                for start in range(0, length, width):
                    end = min(length, start + width)
                    task_list.append([f_bam, f_vcf, c, start, end])                
        else:
            c = chrom
            length = bam.get_reference_length(c)
            for start in range(0, length, width):
                end = min(length, start + width)
                task_list.append([f_bam, f_vcf, c, start, end])   
                
    result_list = []
    pool = None
    if threads > 1:
        pool = mp.Pool(threads)
    for task in task_list:
        if pool:
            r = pool.apply_async(call_snps, tuple(task))
        else:
            r = call_snps(*task)
        result_list.append(r)
    if pool:
        pool.close()
        pool.join()
        result_list = [r.get() for r in result_list]
        
    with open(f_out, "w+") as fw:
        fw.write("##fileformat=VCFv4.1\n")
        for item in sorted(header["SQ"], key=lambda item: item["SN"]):
            fw.write("##contig=<ID=%s,length=%d>\n" % (item["SN"], item["LN"]))
        fw.write("##INFO=<ID=READS,Number=1,Type=String,Description=\"Support reads for base.\"\n")
        fw.write("##INFO=<ID=CELLS,Number=1,Type=String,Description=\"Support cells for base.\"\n")
        fw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        fw.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample)
        
        for r in result_list:
            for snp in r:
                alts = []
                if snp["base1"] != snp["ref_base"]:
                    alts.append(snp["base1"])
                if snp["base2"] != snp["ref_base"]:
                    alts.append(snp["base2"])
                alts = list(sorted(set(alts)))
                alleles = [snp["ref_base"]] + alts
                
                s1 = ",".join(["%s:%d" % (base, reads) for base, reads in sorted(snp["base_reads"].items())])
                s2 = ",".join(["%s:%d" % (base, cells) for base, cells in sorted(snp["base_cells"].items())])
                info = "READS=%s;CELLS=%s" % (s1, s2)
                
                gt1 = alleles.index(snp["base1"])
                gt2 = alleles.index(snp["base2"])
                if gt1 > gt2:
                    gt1, gt2 = gt2, gt1
                gt = "%d/%d" % (gt1, gt2)
                
                line = "\t".join(map(str, [
                    snp["chrom"], 
                    snp["start"] + 1, 
                    snp["id"], 
                    snp["ref_base"], 
                    ",".join(alts), 
                    ".", # QUAL
                    ".", # FILTER
                    info, 
                    "GT", 
                    gt
                ]))
                fw.write(line + "\n")

if __name__ == "__main__":
    main()