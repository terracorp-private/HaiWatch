#!/usr/bin/env python3

import warnings
warnings.simplefilter(action='ignore', category=Warning)
import pandas as pd
import matplotlib.pyplot as plt
import glob

# .idats and corresponding tumorIDs to read
master_file = pd.read_excel('/home/alpha/programs/python_files/datasets/cnv_and_mut/refSet.xlsx')

cnvID =  master_file['txt_idat']
tumorID = master_file['ING_ID']

path_cnv = '/home/alpha/programs/python_files/datasets/cnv_and_mut/k27_cnv/'
path_var = '/home/alpha/programs/python_files/datasets/cnv_and_mut/k27_snp/'



def filterCNV(del_cutoff,gain_cutoff):

    # filter parameters
    ALL_FREQ = 0.05
    READ_DP = 100

    
    # filter CNVP
    cnv = cnvFrame[['chrom','feature','start','end','metrics']]
    cnv = cnv.replace(regex=r'[a-z]',value = '')
    cnv = cnv.replace(regex=r'X',value = 23)
    cnv = cnv.replace(regex=r'Y',value = 24)
    cnvChrom = cnv['chrom']
    cnvStart = cnv['start']
    cnvEnd = cnv['end']
    cnv['cnvCombPos'] = (cnvStart + (cnvEnd - cnvStart) / 2).astype('int')
    cnv['cnvCombPos'] = cnvChrom.astype('str') + ':' + cnv['cnvCombPos'].astype('str')
    cnvFiltered = cnv[(cnv['metrics'] < del_cutoff) | (cnv['metrics'] > gain_cutoff)]

    # filter variants
    vars = varFrame[['Chr','Start','Gene.refGene','ExonicFunc.refGene','Ref','Alt','AAChange.refGene','AF','DP']]
    vars = vars.replace(regex=r'X',value = 23)
    vars = vars.replace(regex=r'Y',value = 24)
    vars = vars[vars['AF'] >= ALL_FREQ]
    vars = vars[vars['DP'] >= READ_DP]
    varChrom = vars['Chr']
    varPos = vars['Start']
    vars['varCombPos'] = varChrom.astype('str') + ':' + varPos.astype('str')
    varIndex = vars['varCombPos']
    
    varRel = vars[vars['ExonicFunc.refGene'] != 'synonymous SNV']

    # iterate through only relevant mutations and whole cnv. Make intersection
    cnvPosList = []
    for index, row in cnvFiltered.iterrows():
        startCNV = row['start']
        endCNV = row['end']
        chromCNV = row['chrom']
        currentDF = varRel[(varRel["Start"].between(startCNV, endCNV)) & (varRel["Chr"] == chromCNV)]
        currentDF['info'] = currentDF['Gene.refGene'].astype('str') + ':' + currentDF['Chr'].astype('str') + ':' + currentDF['Start'].astype('str') + ':' + currentDF['Ref'] + ':' + currentDF['Alt'] + ':' + currentDF['AAChange.refGene'].astype('str')
    
        if len(currentDF) > 0:
            cnvPosList.extend(currentDF['info'].tolist())
    return cnvPosList
    # print(cnvPosList,len(cnvPosList))

# read a .idat and tumorID pair into a dataframe
resultList = []

for f_cnv,f_var in zip(cnvID,tumorID):
    cnv_file = glob.glob(path_cnv + f_cnv + '.bins.igv')[0]
    var_file = glob.glob(path_var + str(f_var) + '-DNA-FFPE_MPILEUP_SNP.recode_filtered_DP15_AF0.1.hg19_multianno._sort.csv')[0]
    cnvFrame = pd.read_csv(cnv_file, sep="\t", header=0, names=['chrom', 'start', 'end', 'feature', 'metrics'])
    varFrame = pd.read_csv(var_file, sep=';')
    resultList.extend(filterCNV(-0.3,0.3))


# print(resultList)

geneNameList = []

for gene in resultList:
    geneName = gene.split(':')[0]
    geneNameList.append(geneName)

df = pd.DataFrame({'geneName' : geneNameList , 'geneTranscript' : resultList})

print(df['geneName'].value_counts())
print(df['geneTranscript'].value_counts())

