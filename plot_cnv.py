#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

# load look up table 
TUMOR = input("gib den Namen der Entität ein: ")
keys = pd.read_excel("/home/alpha/programs/python_files/datasets/cnv_and_mut/refSet" + TUMOR + ".xlsx")
ID = keys["ING_ID"]
CNV_NEEDED = int(input("gib Material-ID ein: "))

# find a .igv file via material id
row_with_cnv = keys[ID == CNV_NEEDED]
cnv_name = row_with_cnv["txt_idat"].to_string(index=False)
cnv_file_path = "/home/alpha/programs/python_files/datasets/cnv_and_mut/" + TUMOR + "_cnv/" + cnv_name + ".bins.igv"
cnv_file = pd.read_csv(cnv_file_path,sep="\t",names=["chr","start","end","feature","value"],header=0)

# standarization
x = cnv_file["feature"]
y = cnv_file["value"].astype("float")


# # hier möchte ich prüfen, ob ich einen Bin < 0 mit auch Nachbarn links und rechts
# # < 0 filtern kann
# helper_neighbor = []
# neighbor = []
# for index,row in cnv_file.iterrows():
#     xVal = row["value"].astype("float")
#     helper_neighbor.append(xVal)


# plot cnv
plt.scatter(x,y,s=1,alpha=0.7)
plt.axhline(y=0)
plt.show()

    



