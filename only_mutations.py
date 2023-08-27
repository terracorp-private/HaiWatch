#!/usr/bin/env python3

import warnings
warnings.simplefilter(action='ignore', category=Warning)
import pandas as pd
import glob

TUMOR = 'mng'

path_cnv = '/home/alpha/programs/python_files/datasets/cnv_and_mut/'+ TUMOR + '_cnv/*'
path_var = '/home/alpha/programs/python_files/datasets/cnv_and_mut/'+ TUMOR + '_snp/*'

MASTER_FILE = pd.read_excel("/home/alpha/programs/python_files/datasets/cnv_and_mut/refSetMNG.back.xlsx")
cnv_files = ";".join(glob.glob(path_cnv))
snv_files = ";".join(glob.glob(path_var))

def sourceIntegrity():
    # check whether the files exsit and create a new dataframe with only existing files
    for index, row in MASTER_FILE.iterrows():
        materialID = str(row["ING_MAT_ID"])
        sentrix = str(row["txt_idat"])
        if sentrix in cnv_files:
            MASTER_FILE.at[index, "cnv_file_exsists"] = "yes"
        else:
            MASTER_FILE.at[index, "cnv_file_exsists"] = "no"

        if materialID in snv_files:
            MASTER_FILE.at[index, "ngs_file_exists"] = "yes"
        else:
            MASTER_FILE.at[index, "ngs_file_exists"] = "no"

    cases = MASTER_FILE[MASTER_FILE["ngs_file_exists"] == "yes"]
    cases = cases[cases["cnv_file_exsists"] == "yes"]

    cnvID =  cases['txt_idat']
    tumorID = cases['ING_MAT_ID']

    return cnvID,tumorID

def filterCNV(cnvFrame,varFrame,ALL_FREQ=0.1,READ_DP=15,del_cutoff=-0.3,gain_cutoff=0.3):

    # filter CNVP
    cnv = cnvFrame[['chrom','feature','start','end','metrics']]
    cnv = cnv.replace(regex=r'[a-z]',value = '')
    cnv = cnv.replace(regex=r'X',value = 23)
    cnv = cnv.replace(regex=r'Y',value = 24)
    cnvFlat = cnv[(cnv['metrics'] > del_cutoff) & (cnv['metrics'] < gain_cutoff)]

    # filter variants
    vars = varFrame[['Chr','Start','Func.refGene','Gene.refGene','ExonicFunc.refGene','Ref','Alt','AAChange.refGene','AF','DP']]
    vars = vars[vars['Func.refGene'] == 'exonic'] # include only exonic variants
    vars = vars.replace(regex=r'X',value = 23) # rename chromosome X to 23
    vars = vars.replace(regex=r'Y',value = 24) # rename chromosome Y to 24
    vars = vars[vars['AF'] >= ALL_FREQ]
    vars = vars[vars['DP'] >= READ_DP]
    varRel = vars[vars['ExonicFunc.refGene'] != 'synonymous SNV']

    mutations = []
    for index, row in cnvFlat.iterrows():
        startCNV = row['start']
        endCNV = row['end']
        chromCNV = row['chrom']
        currentDF = varRel[(varRel["Start"] < startCNV) | (varRel["Start"] > endCNV) & (varRel["Chr"] == chromCNV)]
        currentDF['info'] = currentDF['Gene.refGene'].astype('str') + ':' + currentDF['Chr'].astype('str') + ':' + currentDF['Start'].astype('str') + ':' + currentDF['Ref'] + ':' + currentDF['Alt'] # + ':' + currentDF['AAChange.refGene'].astype('str')
    
        if len(currentDF) > 0:
            mutations.extend(currentDF['info'].tolist())
            
    return mutations

def geneNames(df,gene,f_var):

    if len(gene) != 0:
        df_helper = pd.DataFrame()
        df_helper["gene"] = gene
        df_helper["sample"] = f_var
        df_helper = df_helper.drop_duplicates(subset=['gene'])
        df = pd.concat([df,df_helper])
        
    return df

geneAmplList = []
geneLossList = []

dfMut = pd.DataFrame()

# return only existing files
cnvID,tumorID = sourceIntegrity()

# read a .idat and tumorID pair into a dataframe
for f_cnv,f_var in zip(cnvID,tumorID):
    cnv_file = glob.glob(path_cnv + f_cnv + '.bins.igv')[0]
    var_file = glob.glob(path_var + str(f_var) + '-DNA-FFPE_MPILEUP_SNP.recode_filtered_DP15_AF0.1.hg19_multianno._sort.csv')[0]
    cnvFrame = pd.read_csv(cnv_file, sep="\t", header=0, names=['chrom', 'start', 'end', 'feature', 'metrics'])
    varFrame = pd.read_csv(var_file, sep=';')
    onlyMut = filterCNV(cnvFrame,varFrame)
    dfMut = geneNames(dfMut,onlyMut,f_var)
    print('mutation file : ',var_file,' is ready')

# drop unique values
dfMut = dfMut[dfMut.duplicated(subset=['gene'], keep=False)]
print(dfMut)
dfMut.to_csv(TUMOR + "_only_mut_genes.csv")
