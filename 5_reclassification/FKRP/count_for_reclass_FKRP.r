# Rscript --vanilla count_for_reclass_FKRP.r gnomad_fkrp-1_annotated_v4_confidence_v2.tsv > count_FKRP

#Rscript --vanilla reclassification_plot_FKRP.r count_FKRP


library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

data = read.csv(args[1],header=T,sep="\t")


# Inflection point between the ClinVar P/LP and ClinVar B/LB densities
B2M=-0.7436399

# Use real patient data (8 cohorts) to decide
# M2I was determined according to the peak point of the density plot of the Mild cases
# I2S was determined according to the peak point of the density plot of the Severe cases

# the real patient data used biallelic scores
# -1 was used to calculate the equivalent monoallelic score

M2I=-0.7358121-1
I2S=-1.307241-1



PBH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= B2M & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PBM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= B2M & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PBL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= B2M & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PMH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= M2I & log2(functional_score) < B2M & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PMM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= M2I  & log2(functional_score) < B2M & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PML=subset(data, Confidence_class=="LOW" & log2(functional_score) >= M2I & log2(functional_score) < B2M & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PIH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= I2S & log2(functional_score) < M2I & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PIM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= I2S  & log2(functional_score) < M2I & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PIL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= I2S & log2(functional_score) < M2I & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PSH=subset(data, Confidence_class=="HIGH"  & log2(functional_score) < I2S & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PSM=subset(data, Confidence_class=="MEDIUM"   & log2(functional_score) < I2S & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PSL=subset(data, Confidence_class=="LOW"  & log2(functional_score) < I2S & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))

VBH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= B2M & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VBM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= B2M & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VBL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= B2M & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VMH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= M2I & log2(functional_score) < B2M & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VMM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= M2I  & log2(functional_score) < B2M & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VML=subset(data, Confidence_class=="LOW" & log2(functional_score) >= M2I & log2(functional_score) < B2M & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VIH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= I2S & log2(functional_score) < M2I & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VIM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= I2S  & log2(functional_score) < M2I & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VIL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= I2S & log2(functional_score) < M2I & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VSH=subset(data, Confidence_class=="HIGH"  & log2(functional_score) < I2S & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VSM=subset(data, Confidence_class=="MEDIUM"   & log2(functional_score) < I2S & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VSL=subset(data, Confidence_class=="LOW"  & log2(functional_score) < I2S & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))

BBH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= B2M & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BBM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= B2M & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BBL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= B2M & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BMH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= M2I & log2(functional_score) < B2M & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BMM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= M2I  & log2(functional_score) < B2M & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BML=subset(data, Confidence_class=="LOW" & log2(functional_score) >= M2I & log2(functional_score) < B2M & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BIH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= I2S & log2(functional_score) < M2I & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BIM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= I2S  & log2(functional_score) < M2I & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BIL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= I2S & log2(functional_score) < M2I & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BSH=subset(data, Confidence_class=="HIGH"  & log2(functional_score) < I2S & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BSM=subset(data, Confidence_class=="MEDIUM"   & log2(functional_score) < I2S & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BSL=subset(data, Confidence_class=="LOW"  & log2(functional_score) < I2S & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))


UBH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= B2M  & is.na(Clinvar_clinical_significance))
UBM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= B2M  & is.na(Clinvar_clinical_significance))
UBL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= B2M  & is.na(Clinvar_clinical_significance))
UMH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= M2I & log2(functional_score) < B2M  & is.na(Clinvar_clinical_significance))
UMM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= M2I  & log2(functional_score) < B2M  & is.na(Clinvar_clinical_significance))
UML=subset(data, Confidence_class=="LOW" & log2(functional_score) >= M2I & log2(functional_score) < B2M  & is.na(Clinvar_clinical_significance))
UIH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= I2S & log2(functional_score) < M2I  & is.na(Clinvar_clinical_significance))
UIM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= I2S  & log2(functional_score) < M2I  & is.na(Clinvar_clinical_significance))
UIL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= I2S & log2(functional_score) < M2I  & is.na(Clinvar_clinical_significance))
USH=subset(data, Confidence_class=="HIGH"  & log2(functional_score) < I2S  & is.na(Clinvar_clinical_significance))
USM=subset(data, Confidence_class=="MEDIUM"   & log2(functional_score) < I2S  & is.na(Clinvar_clinical_significance))
USL=subset(data, Confidence_class=="LOW"  & log2(functional_score) < I2S  & is.na(Clinvar_clinical_significance))









cat("Clinvar_clinical_significance","\t", "SMuRF_reclassification","\t", "confidence","\t", "Counts", "\n",sep = "")


cat("P/LP","\t", "SMuRF Benign","\t","High","\t", nrow(PBH), "\n",sep = "")
cat("P/LP","\t", "SMuRF Benign","\t","Medium","\t", nrow(PBM), "\n",sep = "")
cat("P/LP","\t", "SMuRF Benign","\t","Low","\t", nrow(PBL), "\n",sep = "")

cat("P/LP","\t", "SMuRF Mild","\t","High","\t", nrow(PMH), "\n",sep = "")
cat("P/LP","\t", "SMuRF Mild","\t","Medium","\t", nrow(PMM), "\n",sep = "")
cat("P/LP","\t", "SMuRF Mild","\t","Low","\t", nrow(PML), "\n",sep = "")

cat("P/LP","\t", "SMuRF Intermediate","\t","High","\t", nrow(PIH), "\n",sep = "")
cat("P/LP","\t", "SMuRF Intermediate","\t","Medium","\t", nrow(PIM), "\n",sep = "")
cat("P/LP","\t", "SMuRF Intermediate","\t","Low","\t", nrow(PIL), "\n",sep = "")

cat("P/LP","\t", "SMuRF Severe","\t","High","\t", nrow(PSH), "\n",sep = "")
cat("P/LP","\t", "SMuRF Severe","\t","Medium","\t", nrow(PSM), "\n",sep = "")
cat("P/LP","\t", "SMuRF Severe","\t","Low","\t", nrow(PSL), "\n",sep = "")


cat("VUS","\t", "SMuRF Benign","\t","High","\t", nrow(VBH), "\n",sep = "")
cat("VUS","\t", "SMuRF Benign","\t","Medium","\t", nrow(VBM), "\n",sep = "")
cat("VUS","\t", "SMuRF Benign","\t","Low","\t", nrow(VBL), "\n",sep = "")

cat("VUS","\t", "SMuRF Mild","\t","High","\t", nrow(VMH), "\n",sep = "")
cat("VUS","\t", "SMuRF Mild","\t","Medium","\t", nrow(VMM), "\n",sep = "")
cat("VUS","\t", "SMuRF Mild","\t","Low","\t", nrow(VML), "\n",sep = "")

cat("VUS","\t", "SMuRF Intermediate","\t","High","\t", nrow(VIH), "\n",sep = "")
cat("VUS","\t", "SMuRF Intermediate","\t","Medium","\t", nrow(VIM), "\n",sep = "")
cat("VUS","\t", "SMuRF Intermediate","\t","Low","\t", nrow(VIL), "\n",sep = "")

cat("VUS","\t", "SMuRF Severe","\t","High","\t", nrow(VSH), "\n",sep = "")
cat("VUS","\t", "SMuRF Severe","\t","Medium","\t", nrow(VSM), "\n",sep = "")
cat("VUS","\t", "SMuRF Severe","\t","Low","\t", nrow(VSL), "\n",sep = "")

cat("B/LB","\t", "SMuRF Benign","\t","High","\t", nrow(BBH), "\n",sep = "")
cat("B/LB","\t", "SMuRF Benign","\t","Medium","\t", nrow(BBM), "\n",sep = "")
cat("B/LB","\t", "SMuRF Benign","\t","Low","\t", nrow(BBL), "\n",sep = "")

cat("B/LB","\t", "SMuRF Mild","\t","High","\t", nrow(BMH), "\n",sep = "")
cat("B/LB","\t", "SMuRF Mild","\t","Medium","\t", nrow(BMM), "\n",sep = "")
cat("B/LB","\t", "SMuRF Mild","\t","Low","\t", nrow(BML), "\n",sep = "")

cat("B/LB","\t", "SMuRF Intermediate","\t","High","\t", nrow(BIH), "\n",sep = "")
cat("B/LB","\t", "SMuRF Intermediate","\t","Medium","\t", nrow(BIM), "\n",sep = "")
cat("B/LB","\t", "SMuRF Intermediate","\t","Low","\t", nrow(BIL), "\n",sep = "")

cat("B/LB","\t", "SMuRF Severe","\t","High","\t", nrow(BSH), "\n",sep = "")
cat("B/LB","\t", "SMuRF Severe","\t","Medium","\t", nrow(BSM), "\n",sep = "")
cat("B/LB","\t", "SMuRF Severe","\t","Low","\t", nrow(BSL), "\n",sep = "")


cat("Unclassified","\t", "SMuRF Benign","\t","High","\t", nrow(UBH), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Benign","\t","Medium","\t", nrow(UBM), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Benign","\t","Low","\t", nrow(UBL), "\n",sep = "")

cat("Unclassified","\t", "SMuRF Mild","\t","High","\t", nrow(UMH), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Mild","\t","Medium","\t", nrow(UMM), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Mild","\t","Low","\t", nrow(UML), "\n",sep = "")

cat("Unclassified","\t", "SMuRF Intermediate","\t","High","\t", nrow(UIH), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Intermediate","\t","Medium","\t", nrow(UIM), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Intermediate","\t","Low","\t", nrow(UIL), "\n",sep = "")

cat("Unclassified","\t", "SMuRF Severe","\t","High","\t", nrow(USH), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Severe","\t","Medium","\t", nrow(USM), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Severe","\t","Low","\t", nrow(USL), "\n",sep = "")
