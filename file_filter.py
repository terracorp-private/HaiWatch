import pandas as pd
import glob

master_file = pd.read_excel("/home/alpha/programs/python_files/datasets/cnv_and_mut/refSetMNG.back.xlsx")
cnv_files_path = ";".join(glob.glob("/home/alpha/programs/python_files/datasets/cnv_and_mut/mng_cnv/*"))
snv_files_path = ";".join(glob.glob("/home/alpha/programs/python_files/datasets/cnv_and_mut/mng_snp/*"))

cnv_file_helper = []
snv_file_helper = []

for index, row in master_file.iterrows():
    materialID = str(row["ING_MAT_ID"])
    sentrix = str(row["txt_idat"])
    if sentrix in cnv_files_path:
        master_file.at[index, "cnv_file_exsists"] = "yes"
    else:
        master_file.at[index, "cnv_file_exsists"] = "no"

    if materialID in snv_files_path:
        master_file.at[index, "ngs_file_exists"] = "yes"
    else:
        master_file.at[index, "ngs_file_exists"] = "no"

cases = master_file[master_file["ngs_file_exists"] == "yes"]
cases = cases[cases["cnv_file_exsists"] == "yes"]
# cases = master_file[(master_file["cnv_file_exsists"] == "yes") and (master_file["ngs_file_exists"] == "yes")]

print(cases)

