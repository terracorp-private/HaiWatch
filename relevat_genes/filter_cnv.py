#!/usr/bin/env python3

import pandas as pd
import re

reference = pd.read_csv("/home/alpha/programs/python_files/datasets/cnv_and_mut/hg19.ncbiRefSeq.csv")
cnv = pd.read_csv("/home/alpha/programs/python_files/datasets/cnv_and_mut/gbm_cnv/10006823069_R02C01.bins.igv")

reference["seqname"] = reference["seqname"].replace(regex=r'X',value = 23)
reference["seqname"] = reference["seqname"].replace(regex=r'Y',value = 24)
reference["seqname"] = reference["seqname"].replace(regex=r'[_]',value = "")
reference = reference[reference["seqname"] != ""]
reference.dropna()
# reference = reference[reference["seqname"].str.contains("chr.[0-9]{1,2}")]
start = reference["start"].astype("int")
end = reference["end"].astype("int")
gene = reference["gene_name"].astype("str")



print(reference)

