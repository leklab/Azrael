#python3 update_logofunc_revise.py
import csv

def rev(x):
    if x=="A":
        return "T"
    if x=="T":
        return "A"
    if x=="G":
        return "C"
    if x=="C":
        return "G"



old = []
with open("LARGE1_DiMSum_annotated_v6.tsv", 'r') as text:
    lines = text.readlines()
    for line in lines:
        line = line.rstrip()
        old.append(line.split("\t"))

update = []

with open("downloadDataLoGoFunc2.txt", 'r') as text:
    lines = text.readlines()
    for line in lines:
        line = line.rstrip()
        update.append(line.split(","))


update_dict1 = {}
update_dict2 = {}

#key=chrom+hg38_site+wt+variant


for i in update[1:]:
    update_dict1[i[0]+i[1]+i[2]+i[3]] = i[5]
    update_dict2[i[0]+i[1]+i[2]+i[3]] = i[7]

out = [old[0]+['logofunc_pred']+['LoGoFunc_GOF']]



for i in old[1:]:
    #variant_key = '\"19\"' + "\"" + i[1] + "\""+ "\""+ i[4]+ "\""+ "\"" + i[5]+ "\""
    variant_key = '\"22\"' + "\"" + i[2] + "\""+ "\""+ rev(i[5])+ "\""+ "\"" + rev(i[6])+ "\""
    if variant_key in update_dict1.keys():
        troublefix1=update_dict1[variant_key].split("\"")
        troublefix2=update_dict2[variant_key].split("\"")
        logofunc_pred = troublefix1[1]
        LoGoFunc_GOF = troublefix2[1]
    else:
        logofunc_pred = "na"
        LoGoFunc_GOF= "na"
    out.append(i+[logofunc_pred]+[LoGoFunc_GOF])

#with open("Supplementary Table 4_updated_v4.tsv", 'w', newline='') as file:
with open("LARGE1_DiMSum_annotated_v7.tsv", 'w', newline='') as file:
    file_tsv = csv.writer(file, delimiter = "\t")
    file_tsv.writerows(out)
