# Rscript --vanilla count_for_reclass_LARGE1.r gnomad_large1_annotated_v4_confidence_v2.tsv > count_LARGE1

#Rscript --vanilla reclassification_plot_LARGE1.r count_LARGE1


library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

data = read.csv(args[1],header=T,sep="\t")


# Inflection point between the ClinVar P/LP and ClinVar B/LB densities
benign_to_pathogenic=-1.405088


PBH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= benign_to_pathogenic & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PBM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= benign_to_pathogenic & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PBL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= benign_to_pathogenic & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PPH=subset(data, Confidence_class=="HIGH" & log2(functional_score) < benign_to_pathogenic & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PPM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) < benign_to_pathogenic & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))
PPL=subset(data, Confidence_class=="LOW" & log2(functional_score) < benign_to_pathogenic & (Clinvar_clinical_significance == "Likely pathogenic" | Clinvar_clinical_significance == "Pathogenic" | Clinvar_clinical_significance == "Pathogenic/Likely pathogenic"))


VBH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= benign_to_pathogenic & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VBM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= benign_to_pathogenic & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VBL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= benign_to_pathogenic & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VPH=subset(data, Confidence_class=="HIGH" & log2(functional_score) < benign_to_pathogenic & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VPM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) < benign_to_pathogenic & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))
VPL=subset(data, Confidence_class=="LOW" & log2(functional_score) < benign_to_pathogenic & (Clinvar_clinical_significance == "Uncertain significance" | Clinvar_clinical_significance == "Conflicting interpretations of pathogenicity" | Clinvar_clinical_significance == "no interpretation for the single variant"))


BBH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= benign_to_pathogenic & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BBM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= benign_to_pathogenic & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BBL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= benign_to_pathogenic & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BPH=subset(data, Confidence_class=="HIGH" & log2(functional_score) < benign_to_pathogenic & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BPM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) < benign_to_pathogenic & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))
BPL=subset(data, Confidence_class=="LOW" & log2(functional_score) < benign_to_pathogenic & (Clinvar_clinical_significance == "Likely benign" | Clinvar_clinical_significance == "Benign" | Clinvar_clinical_significance == "Benign/Likely benign"))

UBH=subset(data, Confidence_class=="HIGH" & log2(functional_score) >= benign_to_pathogenic  & is.na(Clinvar_clinical_significance))
UBM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) >= benign_to_pathogenic  & is.na(Clinvar_clinical_significance))
UBL=subset(data, Confidence_class=="LOW" & log2(functional_score) >= benign_to_pathogenic  & is.na(Clinvar_clinical_significance))
UPH=subset(data, Confidence_class=="HIGH" & log2(functional_score) < benign_to_pathogenic  & is.na(Clinvar_clinical_significance))
UPM=subset(data, Confidence_class=="MEDIUM" & log2(functional_score) < benign_to_pathogenic  & is.na(Clinvar_clinical_significance))
UPL=subset(data, Confidence_class=="LOW" & log2(functional_score) < benign_to_pathogenic  & is.na(Clinvar_clinical_significance))







cat("Clinvar_clinical_significance","\t", "SMuRF_reclassification","\t", "confidence","\t", "Counts", "\n",sep = "")


cat("P/LP","\t", "SMuRF Benign","\t","High","\t", nrow(PBH), "\n",sep = "")
cat("P/LP","\t", "SMuRF Benign","\t","Medium","\t", nrow(PBM), "\n",sep = "")
cat("P/LP","\t", "SMuRF Benign","\t","Low","\t", nrow(PBL), "\n",sep = "")
cat("P/LP","\t", "SMuRF Pathogenic","\t","High","\t", nrow(PPH), "\n",sep = "")
cat("P/LP","\t", "SMuRF Pathogenic","\t","Medium","\t", nrow(PPM), "\n",sep = "")
cat("P/LP","\t", "SMuRF Pathogenic","\t","Low","\t", nrow(PPL), "\n",sep = "")

cat("VUS","\t", "SMuRF Benign","\t","High","\t", nrow(VBH), "\n",sep = "")
cat("VUS","\t", "SMuRF Benign","\t","Medium","\t", nrow(VBM), "\n",sep = "")
cat("VUS","\t", "SMuRF Benign","\t","Low","\t", nrow(VBL), "\n",sep = "")
cat("VUS","\t", "SMuRF Pathogenic","\t","High","\t", nrow(VPH), "\n",sep = "")
cat("VUS","\t", "SMuRF Pathogenic","\t","Medium","\t", nrow(VPM), "\n",sep = "")
cat("VUS","\t", "SMuRF Pathogenic","\t","Low","\t", nrow(VPL), "\n",sep = "")

cat("B/LB","\t", "SMuRF Benign","\t","High","\t", nrow(BBH), "\n",sep = "")
cat("B/LB","\t", "SMuRF Benign","\t","Medium","\t", nrow(BBM), "\n",sep = "")
cat("B/LB","\t", "SMuRF Benign","\t","Low","\t", nrow(BBL), "\n",sep = "")
cat("B/LB","\t", "SMuRF Pathogenic","\t","High","\t", nrow(BPH), "\n",sep = "")
cat("B/LB","\t", "SMuRF Pathogenic","\t","Medium","\t", nrow(BPM), "\n",sep = "")
cat("B/LB","\t", "SMuRF Pathogenic","\t","Low","\t", nrow(BPL), "\n",sep = "")

cat("Unclassified","\t", "SMuRF Benign","\t","High","\t", nrow(UBH), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Benign","\t","Medium","\t", nrow(UBM), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Benign","\t","Low","\t", nrow(UBL), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Pathogenic","\t","High","\t", nrow(UPH), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Pathogenic","\t","Medium","\t", nrow(UPM), "\n",sep = "")
cat("Unclassified","\t", "SMuRF Pathogenic","\t","Low","\t", nrow(UPL), "\n",sep = "")
