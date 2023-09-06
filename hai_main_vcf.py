#!/usr/bin/env python3

import warnings
warnings.simplefilter(action='ignore', category=Warning)
import pandas as pd
import glob

TUMOR = input("name of the entity ")
THRESHOLD = float(input("loss/ampl threshold "))

# path_cnv = '/home/alpha/programs/python_files/datasets/cnv_and_mut/'+ TUMOR + '_cnv/*'
path_cnv = '/home/alpha/programs/python_files/datasets/cnv_and_mut/'+ TUMOR + '_cnv/*'
path_var = '/home/alpha/programs/python_files/datasets/cnv_and_mut/'+ TUMOR + '_snp/*'

MASTER_FILE = pd.read_excel("/home/alpha/programs/python_files/datasets/cnv_and_mut/refSet" + TUMOR + ".xlsx")
cnv_files = ";".join(glob.glob(path_cnv))
snv_files = ";".join(glob.glob(path_var))

def sourceIntegrity():
    # check whether the files exist and create a new dataframe with only existing files
    for index, row in MASTER_FILE.iterrows():
        materialID = str(row["ING_MAT_ID"])
        sentrix = str(row["txt_idat"])
        if sentrix in cnv_files:
            MASTER_FILE.at[index, "cnv_file_exsists"] = "yes"
        else:
            MASTER_FILE.at[index, "cnv_file_exsists"] = "no"
            print(sentrix,"cnv file not found")

        if materialID in snv_files:
            MASTER_FILE.at[index, "ngs_file_exists"] = "yes"
        else:
            MASTER_FILE.at[index, "ngs_file_exists"] = "no"
            print(materialID,"snv file not found")

    cases = MASTER_FILE[MASTER_FILE["ngs_file_exists"] == "yes"]
    cases = cases[cases["cnv_file_exsists"] == "yes"]

    cnvID =  cases['txt_idat']
    tumorID = cases['ING_MAT_ID']

    return cnvID,tumorID


def filterCNV(cnvFrame,varFrame,ALL_FREQ=0.3,READ_DP=15,del_cutoff=-0.1,gain_cutoff=0.3):

    # filter CNVP
    cnv = cnvFrame[['chrom','feature','start','end','metrics']]
    cnv["chrom"] = cnv["chrom"].replace(regex=r'[a-z]',value = '')
    cnv["chrom"] = cnv["chrom"].replace(regex=r'X',value = 23)
    cnv["chrom"] = cnv["chrom"].replace(regex=r'Y',value = 24)
    cnv["chrom"] = cnv["chrom"].astype("int")
    cnvAmpl = cnv[cnv['metrics'].astype("float") > gain_cutoff]
    # cnvLoss = cnv[(cnv['metrics'] < del_cutoff) & (cnv['metrics'] > -0.45)]
    cnvLoss = cnv[cnv['metrics'].astype("float") < del_cutoff]

    # # Debug. is the position of interest filtered out?
    # cnvLossDebug = cnvLoss[cnvLoss["feature"] == "chr22-0106"]
    # print(cnvLossDebug)


    # filter variants
    vars = varFrame[['Chr','Start','Func.refGene','Gene.refGene','ExonicFunc.refGene','Ref','Alt','AAChange.refGene','AF','DP']]
    vars = vars[vars['Func.refGene'] != 'intronic'] # include only exonic variants
    vars["Chr"] = vars["Chr"].replace(regex=r'X',value = 23) # rename chromosome X to 23
    vars["Chr"] = vars["Chr"].replace(regex=r'Y',value = 24) # rename chromosome Y to 24
    vars["Chr"] = vars["Chr"].astype("int")
    vars = vars[vars['AF'].astype("float") >= ALL_FREQ]
    vars = vars[vars['DP'].astype("int") >= READ_DP]
    varRel = vars[vars['ExonicFunc.refGene'] != 'synonymous SNV']

    # # Debug. is the position of interest filtered out?
    # varRelDebug = varRel[varRel["Start"] == 30051625]
    # print(varRelDebug)


    # iterate through only relevant mutations and whole cnv. Make intersection
    cnvAmplList = []
    for index, row in cnvAmpl.iterrows():
        startCNV = row['start']
        endCNV = row['end']
        chromCNV = row['chrom']
        currentDF = varRel[(varRel["Start"].between(startCNV, endCNV)) & (varRel["Chr"] == chromCNV)]
        currentDF['info'] = currentDF['Gene.refGene'].astype('str') + ':' + currentDF['Chr'].astype('str') + ':' + currentDF['Start'].astype('str') + ':' + currentDF['Ref'] + ':' + currentDF['Alt'] # + ':' + currentDF['AAChange.refGene'].astype('str')
    
        if len(currentDF) > 0:
            cnvAmplList.extend(currentDF['info'].tolist())


    cnvLossList = []
    for index, row in cnvLoss.iterrows():
        startCNV = row['start']
        endCNV = row['end']
        chromCNV = row['chrom']
        currentDF = varRel[(varRel["Start"].between(startCNV, endCNV)) & (varRel["Chr"] == chromCNV)]
        currentDF['info'] = currentDF['Gene.refGene'].astype('str') + ':' + currentDF['Chr'].astype('str') + ':' + currentDF['Start'].astype('str') + ':' + currentDF['Ref'] + ':' + currentDF['Alt'] # + ':' + currentDF['AAChange.refGene'].astype('str')
    
        if len(currentDF) > 0:
            cnvLossList.extend(currentDF['info'].tolist())
            
    return cnvAmplList,cnvLossList
    # print(cnvPosList,len(cnvPosList))

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

dfAmpl = pd.DataFrame()
dfLoss = pd.DataFrame()

cnvID,tumorID = sourceIntegrity()

# read a .idat and tumorID pair into a dataframe
for f_cnv,f_var in zip(cnvID,tumorID):
    cnv_file = glob.glob(path_cnv + f_cnv + '.bins.igv')[0]
    var_file = glob.glob(path_var + str(f_var) + '-DNA-FFPE_MPILEUP_SNP.recode_filtered_DP15_AF0.1.hg19_multianno._sort.csv')[0]
    print(var_file)
    cnvFrame = pd.read_csv(cnv_file, sep="\t", header=0, names=['chrom', 'start', 'end', 'feature', 'metrics'])
    varFrame = pd.read_csv(var_file, sep=';')
    geneAmpl, geneLoss = filterCNV(cnvFrame,varFrame,del_cutoff=-THRESHOLD,gain_cutoff=THRESHOLD)
    dfAmpl = geneNames(dfAmpl,geneAmpl,f_var)
    dfLoss = geneNames(dfLoss,geneLoss,f_var)
    print('mutation file : ',var_file,' is ready')


# dfAmpl = dfAmpl[dfAmpl.duplicated(subset=['gene'], keep=False)]
# dfLoss = dfLoss[dfLoss.duplicated(subset=['gene'], keep=False)]
dfLoss.to_csv(TUMOR + "_lost.csv")
dfAmpl.to_csv(TUMOR + "_ampl.csv")
