#last updated 09212023

#run this script using the following commands
# python SMuRF_sco_nor_FKRP.py all_possible_SNVs_with_clinvar_FKRP_04202023 FKRP_high1_variant_counts.tsv FKRP_low1_variant_counts.tsv fkrp_1to1485.fasta high1_wt_blocks.tsv low1_wt_blocks.tsv normalized_FKRP_block1_func 1
# python SMuRF_sco_nor_FKRP.py all_possible_SNVs_with_clinvar_FKRP_04202023 FKRP_high2_variant_counts.tsv FKRP_low2_variant_counts.tsv fkrp_1to1485.fasta high2_wt_blocks.tsv low2_wt_blocks.tsv normalized_FKRP_block2_func 2
# python SMuRF_sco_nor_FKRP.py all_possible_SNVs_with_clinvar_FKRP_04202023 FKRP_high3_variant_counts.tsv FKRP_low3_variant_counts.tsv fkrp_1to1485.fasta high3_wt_blocks.tsv low3_wt_blocks.tsv normalized_FKRP_block3_func 3
# python SMuRF_sco_nor_FKRP.py all_possible_SNVs_with_clinvar_FKRP_04202023 FKRP_high4_variant_counts.tsv FKRP_low4_variant_counts.tsv fkrp_1to1485.fasta high4_wt_blocks.tsv low4_wt_blocks.tsv normalized_FKRP_block4_func 4
# python SMuRF_sco_nor_FKRP.py all_possible_SNVs_with_clinvar_FKRP_04202023 FKRP_high5_variant_counts.tsv FKRP_low5_variant_counts.tsv fkrp_1to1485.fasta high5_wt_blocks.tsv low5_wt_blocks.tsv normalized_FKRP_block5_func 5
# python SMuRF_sco_nor_FKRP.py all_possible_SNVs_with_clinvar_FKRP_04202023 FKRP_high6_variant_counts.tsv FKRP_low6_variant_counts.tsv fkrp_1to1485.fasta high6_wt_blocks.tsv low6_wt_blocks.tsv normalized_FKRP_block6_func 6



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


#


block_size=250




# read the sequence

inlookuplines=[line.rstrip('\n') for line in inlookup]
inseqouthighlines=[line.rstrip('\n') for line in high]
inseqoutlowlines=[line.rstrip('\n') for line in low]
WThigh=[line.rstrip('\n') for line in WT_high]
WTlow=[line.rstrip('\n') for line in WT_low]
seqline=[line.rstrip('\n') for line in seq]


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
                    hg38_site=site+46755450
                    outFile.write("19"+'\t'+str(hg38_site)+'\t'+str(site)+'\t'+str(codon)+'\t'+WT_nt+'\t'+Variant_nt+'\t'+WT_AA+'\t'+Variant_AA+'\t'+classification+'\t'+str(functional_score)+'\t'+Clinvar_clinical_significance+'\t'+str(highreads)+'\t'+str(lowreads)+'\t'+str(totaltl)+'\t'+str(lowt)+'\n')



inlookup.close()
high.close()
low.close()
seq.close()
WT_high.close()
WT_low.close()
outFile.close()
