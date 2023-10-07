#python3 update_synvep.py
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
#with open("Supplementary Table 4_updated_v3.tsv", 'r') as text:
with open("Supplementary Table 5_updated_v3.tsv", 'r') as text:
    lines = text.readlines()
    for line in lines:
        line = line.rstrip()
        old.append(line.split("\t"))

update = []

#with open("synVep v 1.0.csv", 'r') as text:
with open("synVep v 1.0_LARGE1.csv", 'r') as text:
    lines = text.readlines()
    for line in lines:
        line = line.rstrip()
        update.append(line.split(","))


update_dict = {}

#key=chrom+hg38_site+wt+variant


for i in update[1:]:
    update_dict[i[2]+i[13]+i[4]+i[5]] = i[8]



out = [old[0]+['synvep']]



for i in old[1:]:
    #variant_key = '\"19\"' + "\"" + i[1] + "\""+ "\""+ i[4]+ "\""+ "\"" + i[5]+ "\""
    variant_key = '\"22\"' + "\"" + i[1] + "\""+ "\""+ rev(i[4])+ "\""+ "\"" + rev(i[5])+ "\""
    if variant_key in update_dict.keys():
        troublefix=update_dict[variant_key].split("\"")
        synvep = troublefix[1]
    else:
        synvep = "na"
    out.append(i+[synvep])

#with open("Supplementary Table 4_updated_v4.tsv", 'w', newline='') as file:
with open("Supplementary Table 5_updated_v4.tsv", 'w', newline='') as file:
    file_tsv = csv.writer(file, delimiter = "\t")
    file_tsv.writerows(out)
