#python3

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
#with open("Supplementary Table 4_updated_v2.tsv", 'r') as text:
with open("LARGE1_DiMSum_annotated_v2.tsv", 'r') as text:
    lines = text.readlines()
    for line in lines:
        line = line.rstrip()
        old.append(line.split("\t"))

update = []

#with open("FKRPam.tsv", 'r') as text:
with open("LARGE1am.tsv", 'r') as text:
    lines = text.readlines()
    for line in lines:
        line = line.rstrip()
        update.append(line.split("\t"))


update_dict = {}

for i in update[4:]:
    update_dict[i[0]+i[1]+i[2]+i[3]] = i[8]


out = [old[0]+['am_pathogenicity']]



for i in old[1:]:
    #variant_key = 'chr19' + i[1] + i[4] + i[5]
    variant_key = 'chr22' + i[2] + rev(i[5]) + rev(i[6])
    if variant_key in update_dict.keys():
        am_pathogenicity = update_dict[variant_key]
    else:
        am_pathogenicity = "na"
    out.append(i+[am_pathogenicity])

#with open("Supplementary Table 4_updated_v3.tsv", 'w', newline='') as file:
with open("LARGE1_DiMSum_annotated_v3.tsv", 'w', newline='') as file:
    file_tsv = csv.writer(file, delimiter = "\t")
    file_tsv.writerows(out)
