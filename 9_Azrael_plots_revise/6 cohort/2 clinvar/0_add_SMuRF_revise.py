
#python 0_add_SMuRF_revise.py FKRP_DiMSum_annotated.tsv restructured_data report_and_SMuRF_high_confidence_dimsum

#Rscript --vanilla 1_bar.r report_and_SMuRF

#Rscript --vanilla 2_onset.r report_and_SMuRF

#Rscript --vanilla 3_ck.r report_and_SMuRF


import sys

# input the files

in_tsv=open(sys.argv[1],'rt')
in_pmid=open(sys.argv[2],'rt')
outFile=open(sys.argv[3],'wt')


# read the sequence


in_tsv_lines=[line.rstrip('\n') for line in in_tsv]
in_pmid_lines=[line.rstrip('\n') for line in in_pmid]




# make a header
outFile.write("Patient_no"+'\t'+"Sex"+'\t'+"Age"+'\t'+"Allele1"+'\t'+"Allele2"+'\t'+"Conditions"+'\t'+"Onset_m"+'\t'+"CK_IU_L_ave"+'\t'+"fs1"+'\t'+"fs2"+'\t'+"PMID"+'\t'+"note"+'\n')

for i in range(1, len(in_pmid_lines)):
    this_line=in_pmid_lines[i].split("\t")
    fs1="na"
    fs2="na"
    for j in range(1, len(in_tsv_lines)):
        tltsv=in_tsv_lines[j].split('\t')
        if ("c."+tltsv[2]+tltsv[4]+">"+tltsv[5])==this_line[3] and tltsv[13]=="HIGH":
            fs1=tltsv[41]
        if ("c."+tltsv[2]+tltsv[4]+">"+tltsv[5])==this_line[4] and tltsv[13]=="HIGH":
            fs2=tltsv[41]
    outFile.write(this_line[0]+'\t'+this_line[1]+'\t'+this_line[2]+'\t'+this_line[3]+'\t'+this_line[4]+'\t'+this_line[5]+'\t'+this_line[6]+'\t'+this_line[9]+'\t'+fs1+'\t'+fs2+'\t'+this_line[10]+'\t'+this_line[11]+'\n')









in_tsv.close
in_pmid.close
outFile.close
