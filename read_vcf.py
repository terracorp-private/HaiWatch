import pandas as pd
import re

# vcf = pd.read_csv("/home/alpha/programs/python_files/datasets/cnv_and_mut/294234-DNA-FFPE_MPILEUP.vcf",
               # comment="#",sep="\t")

varFrame = pd.read_csv("/home/alpha/programs/python_files/datasets/cnv_and_mut/mng_snp/78512-DNA-FFPE_PLATYPUS_comp_INDEL.recode_filtered_DP15-1.hg19_multianno._sort.csv",sep=";")

# print(df["AF"].str.replace(",",".").astype("float"))


vars = varFrame[['Chr','Start','Func.refGene','Gene.refGene','ExonicFunc.refGene','Ref','Alt','AAChange.refGene','AF','DP']]

vars = vars[vars["Chr"].str.match("\d")] # remove all non digit-like strings
vars["Chr"] = vars["Chr"].astype("int")

print(vars)
