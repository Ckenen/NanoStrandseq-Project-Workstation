#!/usr/bin/env runsnakemake
import pandas as pd
configfile: "config.yaml" 
THREADS = config["THREADS"]
RUNS = config["HG001_RUNS"] + config["B6D2F1_RUNS"]
DAT = pd.read_excel(config["TABLE"])
DAT.index = DAT["Cell"]
DAT = DAT[DAT["Run"].isin(RUNS)]
DAT["RunCell"] = ["%s/%s" % (run, cell) for run, cell in DAT[["Run", "Cell"]].values]
RUN_CELLS = list(DAT["RunCell"])

def get_species(cell):
    return DAT.loc[cell]["Species"]

def get_strain(cell):
    return DAT.loc[cell]["Strain"]

def get_fasta(cell):
    return config["%s_FASTA" % get_species(cell).upper()]

def get_mmi(cell):
    return config["%s_MMI" % get_species(cell).upper()]

def get_snp_vcf(cell):
    return config["%s_SNP_VCF" % get_strain(cell).upper()]
    
def get_snp_bed(cell):
    return config["%s_SNP_BED" % get_strain(cell).upper()]