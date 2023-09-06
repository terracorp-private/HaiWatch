#!/usr/bin/env python3

import warnings
warnings.simplefilter(action='ignore', category=Warning)
import pandas as pd # der import bitte vor die warnings

TUMOR = input("name of the entity ") # würde folgendes vorschlagen: szs.argv[1]
INDEL = input("for INDEL choose 'y' (dafaul no) ")
GENE = ""

if INDEL == "y":
    dfLoss = pd.read_csv("/home/alpha/programs/python_files/datasets/cnv_and_mut/" + TUMOR + "_lost_INDEL.csv")
    dfAmpl = pd.read_csv("/home/alpha/programs/python_files/datasets/cnv_and_mut/" + TUMOR + "_ampl_INDEL.csv")
else:
    dfLoss = pd.read_csv("/home/alpha/programs/python_files/datasets/cnv_and_mut/" + TUMOR + "_lost.csv")
    dfAmpl = pd.read_csv("/home/alpha/programs/python_files/datasets/cnv_and_mut/" + TUMOR + "_ampl.csv")
# dfMut = pd.read_csv("/home/alpha/programs/python_files/cnvp_and_mutations/" + TUMOR + "_only_mut.csv")


lost_genes = dfLoss['gene'].value_counts() # stellst du hiermit nochmal was an?
lost_case = dfLoss['sample'].value_counts()
ampl_genes = dfAmpl['gene'].value_counts()
# mut_genes = dfMut['gene'].value_counts()

dfLoss["gene_name"] = dfLoss["gene"].str.split(":").str[0]
onlyLostNames = dfLoss.drop_duplicates(subset=["gene_name","sample"]).value_counts("gene_name")

dfAmpl["gene_name"] = dfAmpl["gene"].str.split(":").str[0]
onlyAmplNames = dfAmpl.drop_duplicates(subset=["gene_name","sample"]).value_counts("gene_name")

# print("only name of lost gene ",sortedGenes.value_counts().head(50))
print("samples with lost genes genes ", onlyLostNames.head(20))
print("samples with gained genes names",onlyAmplNames.head(20))
print("Lost genes", lost_genes.head(20))
print("Gained genes", ampl_genes.head(20))


