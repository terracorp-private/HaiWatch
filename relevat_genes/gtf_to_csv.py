#!/usr/bin/env python3

import pandas as pd
from gtfparse import read_gtf

reference = read_gtf("/home/alpha/programs/python_files/datasets/cnv_and_mut/hg19.ncbiRefSeq.gtf")

reference = reference[["seqname","start","end","gene_name"]]

reference.to_csv("/home/alpha/programs/python_files/datasets/cnv_and_mut/hg19.ncbiRefSeq.csv")
