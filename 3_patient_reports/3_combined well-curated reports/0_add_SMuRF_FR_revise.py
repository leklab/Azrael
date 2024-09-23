#python 0_add_SMuRF_FR_revise.py FKRP_DiMSum_annotated_NEW.tsv restructured FKRPregstry_report_and_SMuRF_diumsum_NEW


import sys

# input the files

in_tsv=open(sys.argv[1],'rt')
in_pmid=open(sys.argv[2],'rt')
outFile=open(sys.argv[3],'wt')


# read the sequence


in_tsv_lines=[line.rstrip('\n') for line in in_tsv]
in_pmid_lines=[line.rstrip('\n') for line in in_pmid]



# make a header
outFile.write("Patient_no"+'\t'+"Allele1"+'\t'+"Allele2"+'\t'+"Onset_year"+'\t'+"10MWFDS"+'\t'+"NSAD20itemaddition"+'\t'+"fs1"+'\t'+"fs2"+'\n')

for i in range(1, len(in_pmid_lines)):
    this_line=in_pmid_lines[i].split("\t")
    fs1="na"
    fs2="na"
    for j in range(1, len(in_tsv_lines)):
        tltsv=in_tsv_lines[j].split('\t')
        if tltsv[13]=="HIGH" and ("c."+tltsv[2]+tltsv[4]+">"+tltsv[5])==this_line[11]:
            fs1=tltsv[15]
        if tltsv[13]=="HIGH" and ("c."+tltsv[2]+tltsv[4]+">"+tltsv[5])==this_line[12]:
            fs2=tltsv[15]
    outFile.write(this_line[0]+'\t'+this_line[11]+'\t'+this_line[12]+'\t'+this_line[10]+'\t'+this_line[33]+'\t'+this_line[34]+'\t'+fs1+'\t'+fs2+'\n')









in_tsv.close
in_pmid.close
outFile.close
