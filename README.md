**Saturation Mutagenesis-Reinforced Functional Assays (SMuRF)** is an accessible workflow for deep mutational scanning (DMS) studies. The manuscript is currently available on BioRxiv: https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3. 


This workflow involves a 2-way extension PCR-based saturation mutagenesis and high-throughput functional assays that can be evaluated with short-read Next-generation sequencing. We have organized scripts used for SMuRF across three repositories:    

**Balthazar repository** is used to create the oligo library needed for the saturation mutagenesis as well as perform QC for the resulting construct.    
**Gagarmel repository** is used to process the raw data from NGS to quantify the enrichment of the variants.    
**Azrael repository (this repository)** is used to generate and analyze the functional scores, and plot the results.    

**Requirements:**    
Python (tested on v2.7.16)    
R (tested on v4.1.2)    

**R libraries required:**    
ggplot2    
dplyr    
ggsignif    
scales    
RColorBrewer    
pROC    
ggrepel    
reshape2    

*The commands to run the python and R scripts were put to the beginning of the scripts as annotations.*    

The **first job** that Azrael can accomplish is generating the functional scores using the counts generated by Gargamel. The functional score is calculated following:     

Enrichment of a variant (E_var) in a FACS group is calculated as a ratio of the count of the variant (c_var) to the total count (c_total) at the variant site:    
E_var=c_var/c_total    

Enrichment of the WT (E_WT) is calculated separately for each block. E_WT is calculated as a ratio of the number of the reads without variant (r_WT) to the number of the reads with one or no variant (r_clean).    
E_WT=r_WT/r_clean    

Relative enrichment (rE) is a ratio of the enrichment in the high-glycosylation group to the enrichment in the low-glycosylation group:    
rE_var=E_var_high/E_var_low    
rE_WT=E_WT_high/E_WT_low    

The functional score of a variant is calculated as the ratio of its relative enrichment to that of the WT sequence in the corresponding block, and the SMuRF score is calculated as the log2 value of the functional score:    
Functional_score=rE_var/rE_WT    

For plotting, we used SMuRF scores:     
SMuRF=log2(Functional_score)    

**0. Test pipeline improvements:**    
This directory contains the scripts that were used to compare results from different pipelines during development.    

**1. Basic figures:**    
This directory contains the scripts that can be used to generate the following figures:     
1.1 Variant types-SMuRF    
1.2 ClinVar classification-SMuRF    
1.3 Enzymatic domains-SMuRF    

**2. Additional analyses:**    
This directory contains the scripts that can be used to generate the following figures:     
2.1 gnomAD-SMuRF    
2.2 ROC of computational predictors and SMuRF    
2.3 computational predictor score-SMuRF    

**3. Patient reports:**    
This directory contains the scripts that can be used to generate the following figures:     
3.1 Disease severity-SMuRF    
3.2 Onset-SMuRF    
3.3 CK values-SMuRF    

**4. Potential novel discoveries:**    
This directory contains the scripts that can be used to generate the following figures:     
4.1 SynVep scores-SMuRF    
4.2 Phylop-SMuRF    

**5. Reclassification:**     
This directory contains the scripts that can be used to classify and/or reclassify variants based on the SMuRF functional scores and confidence scores.    

**6. ppVSV enrichment:**     
ppVSV enrichment was evaluated with nanopore sequencing provided by Plasmidusaurus. This directory contains scripts to calculate the functional score based on the nanopore sequencing results.   

**7. updates:**    
This directory contains the scripts that can be used for following jobs:    
7.1 Adding scores from a new computational tool to the spreadsheet.    
7.2 LoGoFunc-SMuRF    
7.3 Heatmap for comparing predictors    
7.4 Enzymatic domains-computational scores & ClinVar-computational scores    
7.5 gnomAD was updated to v4; we generated the gnomADv4-SMuRF figures.    
