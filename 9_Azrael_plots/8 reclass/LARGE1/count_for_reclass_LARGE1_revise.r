# Rscript --vanilla count_for_reclass_LARGE1_revise.r LARGE1_DiMSum_annotated_v6.tsv > count_LARGE1_revise



library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

data = read.csv(args[1],header=T,sep="\t")
data = data[data$confidence == "HIGH",]

benign_to_pathogenic=-1.322896



PBH=subset(data,  fitness_normalized >= benign_to_pathogenic & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PPH=subset(data,  fitness_normalized < benign_to_pathogenic & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))


VBH=subset(data,  fitness_normalized >= benign_to_pathogenic & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VPH=subset(data,  fitness_normalized < benign_to_pathogenic & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))


BBH=subset(data,  fitness_normalized >= benign_to_pathogenic & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BPH=subset(data,  fitness_normalized < benign_to_pathogenic & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))

UBH=subset(data,  fitness_normalized >= benign_to_pathogenic  & is.na(Clinvar_clinical_significance))
UPH=subset(data,  fitness_normalized < benign_to_pathogenic  & is.na(Clinvar_clinical_significance))





cat("Clinvar_clinical_significance","\t", "SMuRF_reclassification","\t", "confidence","\t", "Counts", "\n",sep = "")


cat("P/LP","\t", "SMuRF Benign","\t","High","\t", nrow(PBH), "\n",sep = "")
cat("P/LP","\t", "SMuRF Pathogenic","\t","High","\t", nrow(PPH), "\n",sep = "")

cat("VUS","\t", "SMuRF Benign","\t","High","\t", nrow(VBH), "\n",sep = "")
cat("VUS","\t", "SMuRF Pathogenic","\t","High","\t", nrow(VPH), "\n",sep = "")

cat("B/LB","\t", "SMuRF Benign","\t","High","\t", nrow(BBH), "\n",sep = "")
cat("B/LB","\t", "SMuRF Pathogenic","\t","High","\t", nrow(BPH), "\n",sep = "")

cat("Unclassified","\t", "SMuRF Benign","\t","High","\t", nrow(UBH), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Pathogenic","\t","High","\t", nrow(UPH), "\n",sep = "")
