# Rscript --vanilla count_for_reclass_FKRP_revise.r FKRP_DiMSum_annotated_v2.tsv > count_FKRP_revise



library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

data = read.csv(args[1],header=T,sep="\t")
data = data[data$confidence == "HIGH",]

# Inflection point between the ClinVar P/LP and ClinVar B/LB densities
B2M=-0.8727984

# Use real patient data (8 cohorts) to decide
# M2I was determined according to the peak point of the density plot of the Mild cases
# I2S was determined according to the peak point of the density plot of the Severe cases

# the real patient data used biallelic scores
# -1 was used to calculate the equivalent monoallelic score

M2I=-0.6183953-1
I2S=-1.244618-1



PBH=subset(data,  fitness_normalized >= B2M & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PMH=subset(data,  fitness_normalized >= M2I & fitness_normalized < B2M & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PIH=subset(data,  fitness_normalized >= I2S & fitness_normalized < M2I & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PSH=subset(data,  fitness_normalized < I2S & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))

VBH=subset(data,  fitness_normalized >= B2M & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VMH=subset(data,  fitness_normalized >= M2I & fitness_normalized < B2M & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VIH=subset(data,  fitness_normalized >= I2S & fitness_normalized < M2I & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VSH=subset(data, fitness_normalized < I2S & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))

BBH=subset(data,  fitness_normalized >= B2M & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BMH=subset(data,  fitness_normalized >= M2I & fitness_normalized < B2M & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BIH=subset(data,  fitness_normalized >= I2S & fitness_normalized < M2I & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BSH=subset(data,  fitness_normalized < I2S & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))


UBH=subset(data,  fitness_normalized >= B2M  & is.na(Clinvar_clinical_significance))
UMH=subset(data,  fitness_normalized >= M2I & fitness_normalized < B2M  & is.na(Clinvar_clinical_significance))
UIH=subset(data,  fitness_normalized >= I2S & fitness_normalized < M2I  & is.na(Clinvar_clinical_significance))
USH=subset(data,  fitness_normalized < I2S  & is.na(Clinvar_clinical_significance))





cat("Clinvar_clinical_significance","\t", "SMuRF_reclassification","\t", "confidence","\t", "Counts", "\n",sep = "")


cat("P/LP","\t", "SMuRF Benign","\t","High","\t", nrow(PBH), "\n",sep = "")

cat("P/LP","\t", "SMuRF Mild","\t","High","\t", nrow(PMH), "\n",sep = "")

cat("P/LP","\t", "SMuRF Intermediate","\t","High","\t", nrow(PIH), "\n",sep = "")

cat("P/LP","\t", "SMuRF Severe","\t","High","\t", nrow(PSH), "\n",sep = "")

cat("VUS","\t", "SMuRF Benign","\t","High","\t", nrow(VBH), "\n",sep = "")

cat("VUS","\t", "SMuRF Mild","\t","High","\t", nrow(VMH), "\n",sep = "")

cat("VUS","\t", "SMuRF Intermediate","\t","High","\t", nrow(VIH), "\n",sep = "")

cat("VUS","\t", "SMuRF Severe","\t","High","\t", nrow(VSH), "\n",sep = "")

cat("B/LB","\t", "SMuRF Benign","\t","High","\t", nrow(BBH), "\n",sep = "")

cat("B/LB","\t", "SMuRF Mild","\t","High","\t", nrow(BMH), "\n",sep = "")

cat("B/LB","\t", "SMuRF Intermediate","\t","High","\t", nrow(BIH), "\n",sep = "")

cat("B/LB","\t", "SMuRF Severe","\t","High","\t", nrow(BSH), "\n",sep = "")

cat("Unclassified","\t", "SMuRF Benign","\t","High","\t", nrow(UBH), "\n",sep = "")

cat("Unclassified","\t", "SMuRF Mild","\t","High","\t", nrow(UMH), "\n",sep = "")

cat("Unclassified","\t", "SMuRF Intermediate","\t","High","\t", nrow(UIH), "\n",sep = "")

cat("Unclassified","\t", "SMuRF Severe","\t","High","\t", nrow(USH), "\n",sep = "")
