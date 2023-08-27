
import pandas as pd
import glob

TUMOR = input("Which entity do you want to standard? ")
PATH_CNV = glob.glob("/home/alpha/programs/python_files/datasets/cnv_and_mut/" + TUMOR + "_cnv/*")

for file in PATH_CNV:
    cnv_file = pd.read_csv(file,sep="\t",names=["chr","start","end","feature","metrics"],header=0)
    y = cnv_file["metrics"].astype("float")
    yStand = (y - y.mean()) / y.std()
    cnv_file["metrics"] = yStand
    cnv_file.to_csv("/home/alpha/programs/python_files/datasets/cnv_and_mut/mng_cnv_stand/" + file.split("/")[-1],sep="\t",index=False)

