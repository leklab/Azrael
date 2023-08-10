#python clinvar_patient_reports_LARGE1.py gnomad_large1_annotated_v2.tsv clinvar_result_04202023_LARGE1.txt LARGE1_conditions


clip1=7
clip2=89


import sys

# input the files

in_tsv=open(sys.argv[1],'rt')
in_clinvar=open(sys.argv[2],'rt')
outFile=open(sys.argv[3],'wt')


# read the sequence


in_tsv_lines=[line.rstrip('\n') for line in in_tsv]
in_clinvar_lines=[line.rstrip('\n') for line in in_clinvar]




# make a header
outFile.write("chr"+'\t'+"hg38_site"+'\t'+"site"+'\t'+"codon"+'\t'+"WT_nt"+'\t'+"Variant_nt"+'\t'+"WT_AA"+'\t'+"Variant_AA"+'\t'+"classification"+'\t'+"functional_score"+'\t'+"Clinvar_clinical_significance"+'\t'+"Clinvar_conditions"+'\n')





for i in range(1, len(in_tsv_lines)):
    tsv_this_line=in_tsv_lines[i].split("\t")
    chr=tsv_this_line[0]
    hg38_site=tsv_this_line[1]
    site=tsv_this_line[2]
    codon=tsv_this_line[3]
    WT_nt=tsv_this_line[4]
    Variant_nt=tsv_this_line[5]
    WT_AA=tsv_this_line[6]
    Variant_AA=tsv_this_line[7]
    classification=tsv_this_line[8]
    functional_score=tsv_this_line[9]
    Clinvar_clinical_significance=tsv_this_line[10]

    if ((Clinvar_clinical_significance=="Pathogenic") | (Clinvar_clinical_significance=="Pathogenic/likely pathogenic") | (Clinvar_clinical_significance=="Likely pathogenic")):
        for LLL in range(clip1,len(in_clinvar_lines)-clip2):
            cvtl=in_clinvar_lines[LLL].split("\t")
            if len(cvtl)>=5:
                cvnote=cvtl[4].split("(")
                cvtldot=cvtl[0].split("c.")
                if len(cvtldot)>=2:
                    cvtlspw=cvtldot[1].split(">")
                    if len(cvtlspw)>=2:
                        cvtlvnt=cvtlspw[1].split(" (")
                        if cvtlspw[0]==site+WT_nt and cvtlvnt[0]==Variant_nt:
                            if cvtl[3] != "not provided":
                                Clinvar_conditions=cvtl[3]
        multi=Clinvar_conditions.split("|")
        if len(multi)<2:
            outFile.write(chr+'\t'+hg38_site+'\t'+site+'\t'+codon+'\t'+WT_nt+'\t'+Variant_nt+'\t'+WT_AA+'\t'+Variant_AA+'\t'+classification+'\t'+functional_score+'\t'+Clinvar_clinical_significance+'\t'+Clinvar_conditions+'\n')
        else:
            for mc in range(0,len(multi)):
                if multi[mc]!="not provided":
                    outFile.write(chr+'\t'+hg38_site+'\t'+site+'\t'+codon+'\t'+WT_nt+'\t'+Variant_nt+'\t'+WT_AA+'\t'+Variant_AA+'\t'+classification+'\t'+functional_score+'\t'+Clinvar_clinical_significance+'\t'+multi[mc]+'\n')






in_tsv.close
in_clinvar.close
outFile.close
