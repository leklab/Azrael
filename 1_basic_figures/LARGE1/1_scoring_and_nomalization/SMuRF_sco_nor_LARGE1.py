#last updated:09212023

#run this script using the following commands
#python SMuRF_sco_nor_LARGE1.py all_possible_SNVs_with_clinvar_LARGE1_04202023 LARGE1_high1_variant_counts.tsv LARGE1_low1_variant_counts.tsv large1_1to2268.fasta high1_wt_blocks.tsv low1_wt_blocks.tsv normalized_LARGE1_block1_func 1 LARGE1_GRCh38_mapping.tsv
#python SMuRF_sco_nor_LARGE1.py all_possible_SNVs_with_clinvar_LARGE1_04202023 LARGE1_high2_variant_counts.tsv LARGE1_low2_variant_counts.tsv large1_1to2268.fasta high2_wt_blocks.tsv low2_wt_blocks.tsv normalized_LARGE1_block2_func 2 LARGE1_GRCh38_mapping.tsv
#python SMuRF_sco_nor_LARGE1.py all_possible_SNVs_with_clinvar_LARGE1_04202023 LARGE1_high3_variant_counts.tsv LARGE1_low3_variant_counts.tsv large1_1to2268.fasta high3_wt_blocks.tsv low3_wt_blocks.tsv normalized_LARGE1_block3_func 3 LARGE1_GRCh38_mapping.tsv
#python SMuRF_sco_nor_LARGE1.py all_possible_SNVs_with_clinvar_LARGE1_04202023 LARGE1_high4_variant_counts.tsv LARGE1_low4_variant_counts.tsv large1_1to2268.fasta high4_wt_blocks.tsv low4_wt_blocks.tsv normalized_LARGE1_block4_func 4 LARGE1_GRCh38_mapping.tsv
#python SMuRF_sco_nor_LARGE1.py all_possible_SNVs_with_clinvar_LARGE1_04202023 LARGE1_high5_variant_counts.tsv LARGE1_low5_variant_counts.tsv large1_1to2268.fasta high5_wt_blocks.tsv low5_wt_blocks.tsv normalized_LARGE1_block5_func 5 LARGE1_GRCh38_mapping.tsv
#python SMuRF_sco_nor_LARGE1.py all_possible_SNVs_with_clinvar_LARGE1_04202023 LARGE1_high6_variant_counts.tsv LARGE1_low6_variant_counts.tsv large1_1to2268.fasta high6_wt_blocks.tsv low6_wt_blocks.tsv normalized_LARGE1_block6_func 6 LARGE1_GRCh38_mapping.tsv
#python SMuRF_sco_nor_LARGE1.py all_possible_SNVs_with_clinvar_LARGE1_04202023 LARGE1_high7_variant_counts.tsv LARGE1_low7_variant_counts.tsv large1_1to2268.fasta high7_wt_blocks.tsv low7_wt_blocks.tsv normalized_LARGE1_block7_func 7 LARGE1_GRCh38_mapping.tsv
#python SMuRF_sco_nor_LARGE1.py all_possible_SNVs_with_clinvar_LARGE1_04202023 LARGE1_high8_variant_counts.tsv LARGE1_low8_variant_counts.tsv large1_1to2268.fasta high8_wt_blocks.tsv low8_wt_blocks.tsv normalized_LARGE1_block8_func 8 LARGE1_GRCh38_mapping.tsv
#python SMuRF_sco_nor_LARGE1.py all_possible_SNVs_with_clinvar_LARGE1_04202023 LARGE1_high9_variant_counts.tsv LARGE1_low9_variant_counts.tsv large1_1to2268.fasta high9_wt_blocks.tsv low9_wt_blocks.tsv normalized_LARGE1_block9_func 9 LARGE1_GRCh38_mapping.tsv
#python SMuRF_sco_nor_LARGE1.py all_possible_SNVs_with_clinvar_LARGE1_04202023 LARGE1_high10_variant_counts.tsv LARGE1_low10_variant_counts.tsv large1_1to2268.fasta high10_wt_blocks.tsv low10_wt_blocks.tsv normalized_LARGE1_block10_func 10 LARGE1_GRCh38_mapping.tsv



import sys


# input the files

inlookup=open(sys.argv[1],'rt')
high=open(sys.argv[2],'rt')
low=open(sys.argv[3],'rt')
seq=open(sys.argv[4],'rt')
WT_high=open(sys.argv[5],'rt')
WT_low=open(sys.argv[6],'rt')
outFile=open(sys.argv[7],'wt')
block=int(sys.argv[8])
hg38sitefile=open(sys.argv[9],'rt')


#


block_size=227




# read the sequence

inlookuplines=[line.rstrip('\n') for line in inlookup]
inseqouthighlines=[line.rstrip('\n') for line in high]
inseqoutlowlines=[line.rstrip('\n') for line in low]
WThigh=[line.rstrip('\n') for line in WT_high]
WTlow=[line.rstrip('\n') for line in WT_low]
seqline=[line.rstrip('\n') for line in seq]
hg38sitefileline=[line.rstrip('\n') for line in hg38sitefile]


# make a header
outFile.write("chr"+'\t'+"hg38_site"+'\t'+"site"+'\t'+"codon"+'\t'+"WT_nt"+'\t'+"Variant_nt"+'\t'+"WT_AA"+'\t'+"Variant_AA"+'\t'+"classification"+'\t'+"functional_score"+'\t'+"Clinvar_clinical_significance"+'\t'+"high_reads"+'\t'+"low_reads"+'\t'+"high_site_all_reads"+'\t'+"low_site_all_reads"+'\n')


#calculate the WT enrichment
count_in_high=WThigh[0].split("\t")
count_in_low=WTlow[0].split("\t")

WT_in_high=float(count_in_high[1])
Variant_in_high=float(count_in_high[2])

WT_in_low=float(count_in_low[1])
Variant_in_low=float(count_in_low[2])

WT_enr=(WT_in_high/(WT_in_high+Variant_in_high))/(WT_in_low/(WT_in_low+Variant_in_low))


#calculate the high enrichment and the low enrichment
# functional score is high/low

for i in range (1,len(inseqouthighlines)):
    thisline=inseqouthighlines[i].split("\t")
    thislinelow=inseqoutlowlines[i].split("\t")

    Atl=thisline[2]
    Ctl=thisline[3]
    Gtl=thisline[4]
    Ttl=thisline[5]

    lowAtl=thislinelow[2]
    lowCtl=thislinelow[3]
    lowGtl=thislinelow[4]
    lowTtl=thislinelow[5]

    totaltl=float(Atl)+float(Ctl)+float(Gtl)+float(Ttl)
    lowt=float(lowAtl)+float(lowCtl)+float(lowGtl)+float(lowTtl)

    Aenr=float(Atl)/float(totaltl)
    Cenr=float(Ctl)/float(totaltl)
    Genr=float(Gtl)/float(totaltl)
    Tenr=float(Ttl)/float(totaltl)

    Aenrlow=float(lowAtl)/float(lowt)
    Cenrlow=float(lowCtl)/float(lowt)
    Genrlow=float(lowGtl)/float(lowt)
    Tenrlow=float(lowTtl)/float(lowt)

    dict={

      "A": Aenr,
      "T": Tenr,
      "G": Genr,
      "C": Cenr
    }


    dictlow={

      "A": Aenrlow,
      "T": Tenrlow,
      "G": Genrlow,
      "C": Cenrlow
    }

    dicthighreads={

      "A": Atl,
      "T": Ttl,
      "G": Gtl,
      "C": Ctl
    }

    dictlowreads={

      "A": lowAtl,
      "T": lowTtl,
      "G": lowGtl,
      "C": lowCtl
    }





    site=i+(block-1)*block_size

    for j in range((i-1+block_size*(block-1))*3+1,(i-1+block_size*(block-1))*3+4):
        lookuptl=inlookuplines[j].split("\t")
        codon=lookuptl[1]
        WT_nt=lookuptl[2]
        Variant_nt=lookuptl[3]
        Variant_AA="NA"
        if codon!="1":
            Variant_AA=lookuptl[5]
        WT_AA=lookuptl[4]
        classification=lookuptl[6]
        Clinvar_clinical_significance=lookuptl[7]



        for x in ("A","T","G","C"):
            if Variant_nt==x:
                highenrich=dict[x]
                lowenrich=dictlow[x]
                highreads=dicthighreads[x]
                lowreads=dictlowreads[x]

                if lowenrich != 0:
                    functional_score = (highenrich/lowenrich)/WT_enr
                    for mapping in range(1, len(hg38sitefileline)):
                        hg38sitefilethisline=hg38sitefileline[mapping].split("\t")
                        tp_CDS=hg38sitefilethisline[0]
                        tp_hg38_l=hg38sitefilethisline[1].split(":")
                        tp_hg38=tp_hg38_l[1]
                        if str(site)==tp_CDS:
                            hg38_site=tp_hg38
                    outFile.write("22"+'\t'+str(hg38_site)+'\t'+str(site)+'\t'+str(codon)+'\t'+WT_nt+'\t'+Variant_nt+'\t'+WT_AA+'\t'+Variant_AA+'\t'+classification+'\t'+str(functional_score)+'\t'+Clinvar_clinical_significance+'\t'+str(highreads)+'\t'+str(lowreads)+'\t'+str(totaltl)+'\t'+str(lowt)+'\n')



inlookup.close()
high.close()
low.close()
seq.close()
WT_high.close()
WT_low.close()
outFile.close()
hg38sitefile.close()
