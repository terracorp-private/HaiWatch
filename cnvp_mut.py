#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

cnvp_file = pd.read_csv('~/programs/python_files/datasets/cnv_and_mut/207343240093_R05C01.bins.igv',
                        sep="\t",header=0,names=['chrom', 'start', 'end', 'feature', 'metrics'])

variants_file = pd.read_csv('~/programs/python_files/datasets/cnv_and_mut/341874-DNA-FFPE_MPILEUP_SNP.recode_filtered_DP15_AF0.1.hg19_multianno.csv',sep=',')
print(variants_file.head())
print(variants_file.columns)



# filter CNVP
cnv = cnvp_file[['chrom','feature','start','end','metrics']]
cnv = cnv.replace(regex=r'[a-z]',value = '')
cnv = cnv.replace(regex=r'X',value = 23)
cnv = cnv.replace(regex=r'Y',value = 24)
cnvChrom = cnv['chrom']
cnvStart = cnv['start']
cnvEnd = cnv['end']
cnv['cnvCombPos'] = (cnvStart + (cnvEnd - cnvStart) / 2).astype('int')
cnv['cnvCombPos'] = cnvChrom.astype('str') + ':' + cnv['cnvCombPos'].astype('str')
cnvFiltered = cnv[(cnv['metrics'] < -0.3) | (cnv['metrics'] > 0.3)]

print(cnv.head())


# filter variants
vars = variants_file[['Chr','Start','Gene.refGene','ExonicFunc.refGene']]
vars = vars.replace(regex=r'X',value = 23)
vars = vars.replace(regex=r'Y',value = 24)
varChrom = vars['Chr']
varPos = vars['Start']
vars['varCombPos'] = varChrom.astype('str') + ':' + varPos.astype('str')
varIndex = vars['varCombPos']


cnvPosList = []
for index, row in cnvFiltered.iterrows():
    startCNV = row['start']
    endCNV = row['end']
    chromCNV = row['chrom']
    for index, row in vars.iterrows():
        posVAR = row['Start']
        chromVAR = row['Chr']
        gen = row['Gene.refGene']
        func = row['ExonicFunc.refGene']
        if (func == 'nonsynonymous SNV') | (func == 'stopgain'):
            if (startCNV <= posVAR <= endCNV) and (chromVAR == chromCNV):
                #cnvMutInfo = str(chromVAR) + ':' + str(posVAR) + ':' + str(gen) + ':' + func
                cnvMutIntersect = str(chromVAR) + ':' + str(posVAR)
                cnvPosList.append(cnvMutIntersect)
                print(cnvPosList,len(cnvPosList))


yPos = -0.3
mutYPos = []
mutXPos = cnvPosList
for i in mutXPos:
    mutYPos.append(yPos)
    print(mutYPos)
    



#cnvPosList = []
#for index, row in vars.iterrows():
#    posVAR = row['position']
#    chromVAR = row['chrom']
#    for index, row in cnv.iterrows():
#        startCNV = row['start']
#        endCNV = row['end']
#        chromCNV = row['chrom']
#        if (startCNV <= posVAR <= endCNV) and (chromVAR == chromCNV):
#            cnvMutIntersect = str(chromVAR) + ':' + str(posVAR)
#            cnvPosList.append(cnvMutIntersect)
#print(cnvPosList,len(cnvPosList))




## Plot CNVP or vars
plt.title("CNV and SNV")
y_cnv = cnv['metrics']
x_cnv = cnv['feature'].astype('str')
plt.xlabel('thingy which looks like X')
plt.ylabel('some random points up and down')
plt.scatter(x_cnv,y_cnv,s=4,alpha=0.7)
plt.scatter(mutXPos,mutYPos,s=15,alpha=0.8)
plt.axhline(y = 0, color ='k')
plt.show()

## Plot CNVP or vars
#y_cnv = cnv['metrics']
#x_cnv = cnv['chrom']
#plt.scatter(x_cnv,y_cnv,s=0.4,alpha=0.7)
#plt.show()

#y_var = vars['position']
#x_var = vars['chrom']

#plt.scatter(x_var,y_var,s=10,alpha=0.7)
#plt.show()








