**Saturation Mutagenesis-Reinforced Functional Assays (SMuRF)** is an accessible workflow for deep mutational scanning (DMS) studies. The manuscript is currently available on [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3). This workflow involves a 2-way extension PCR-based saturation mutagenesis and high-throughput functional assays that can be evaluated with short-read Next-generation sequencing.

We have organized scripts used for SMuRF across three repositories:    
* [**Balthazar repository**](https://github.com/leklab/Balthazar) is used to create the oligo library needed for the saturation mutagenesis as well as perform QC for the resulting construct.    
* [**Gagarmel repository**](https://github.com/leklab/Gargamel) is used to process the raw data from NGS to quantify the enrichment of the variants.    
* **Azrael repository (this repository)** is used to generate and analyze the functional scores, and plot the results.    

## Requirements
* Python (tested on v2.7.16)    
* R (tested on v4.1.2)    

### R libraries required    
```
ggplot2    
dplyr    
ggsignif    
scales    
RColorBrewer    
pROC    
ggrepel    
reshape2    
```
*The commands to run the python and R scripts were put to the beginning of the scripts as annotations.*    

## SMuRF Functional Score Calculation
The **first job** that Azrael can accomplish is generating the functional scores using the counts generated by Gargamel. The functional score is calculated following:     

Enrichment of a variant (`E_var`) in a FACS group is calculated as a ratio of the count of the variant (`c_var`) to the total count (`c_total`) at the variant site:    
`E_var=c_var/c_total`    

Enrichment of the WT (`E_WT`) is calculated separately for each block. `E_WT` is calculated as a ratio of the number of the reads without variant (`r_WT`) to the number of the reads with one or no variant (`r_clean`).    
`E_WT=r_WT/r_clean`    

Relative enrichment (`rE`) is a ratio of the enrichment in the high-glycosylation group to the enrichment in the low-glycosylation group:    
`rE_var=E_var_high/E_var_low`    
`rE_WT=E_WT_high/E_WT_low`    

The functional score of a variant (in one biological replicate) is calculated as the ratio of its relative enrichment to that of the WT sequence in the corresponding block, and the SMuRF score is calculated as the log2 value of the functional score:    
`Functional_score=rE_var/rE_WT`    

SMuRF score is generated by combing 3 biological replicates (see below).

## Fitness score combination using DiMSum
The DiMSum pipeline was employed to process the variant counts to generate a fitness score across all biological replicates.

For data preprocessing, the input files were prepared by block to account for block-specific WT variant counts. Variants with a count of 0 were replaced by the WT value. The variant counts for each replicate within each gene block were processed through [DiMSum](https://github.com/lehner-lab/DiMSum/) from the STEAM stage. Default settings were used for the pipeline, aside from --sequenceType noncoding to calculate the scores on a variant-level instead of AA-level. The pipeline aggregated these counts to generate a fitness score for each variant. DiMSum incorporates both the frequency of each variant and the observed effects on fitness across all replicates and blocks, producing a merged fitness score. To ensure accurate and reliable fitness estimates, DiMSum employs an advanced error modeling approach. The model addresses multiple sources of error to produce robust fitness scores. Firstly, Poisson-based errors are used to account for the high variance in fitness estimates that arise from low sequencing depth. Since count data often exhibit over-dispersion relative to a Poisson distribution, DiMSum integrates both multiplicative and additive error terms to model this variance effectively. Multiplicative errors, which scale with the variant sequencing counts, address variability arising from workflow-related inconsistencies. In contrast, additive errors are independent of sequencing counts and affect all variants uniformly, accounting for differential handling across replicates. These error components are summarized into a variable, sigma (σ), which is used to adjust fitness estimates on a variant-level. The primary output from the DiMSum pipeline is the 1) “DiMSum score” – fitness score from the direct output of the DiMSum pipeline, converted to the log2 scale, and 2) sigma – a numerical value representing the reliability of the variant count based on error modeling.

## Fitness score normalization
To normalize the fitness scores, we used all synonymous variants as the reference, assuming these variants should have no impact and thus a fitness score of 1. We first transformed the fitness scores of the synonymous variants from their natural logarithm form to a linear scale using the exponential function. For each block, we calculated the mean fitness score of the synonymous variants, treating this mean as the scaling factor. Next, we normalized the fitness scores of all variants by transforming their scores back to the linear scale, dividing by the block's scaling factor, and then applying a log2 transformation to the normalized values to get the “SMuRF score”. This approach ensures that the scores are adjusted relative to the neutral impact of synonymous variants, allowing for consistent comparisons across different blocks.

## Confidence score generation and classification
Following the DiMSum pipeline for variant score combination and error modeling, each identified variant was assigned a sigma (σ) value representing the estimated error. A higher sigma indicates a higher degree of error in the variant assessment. To classify variants into high and low confidence groups, we applied a threshold based on the 1.5 * interquartile range (IQR) rule, a commonly used outlier detection method. First, we calculated the first quartile (Q1) and third quartile (Q3) of the sigma distribution across all variants. The interquartile range (IQR) was then determined as IQR = Q3 - Q1. The upper threshold for high confidence variants was set at Q3 + 1.5*IQR. Any variant with a sigma value exceeding this threshold was classified as low confidence, while variants below the threshold were designated as high confidence. For all downstream analyses, low confidence variants were excluded to ensure the robustness of results. This approach provided an objective method to filter out potentially unreliable sites prior to further examination and interpretation of the data.


## Function score analysis and Figures
**0. Test pipeline improvements:**    
This directory contains the scripts that were used to compare results from different pipelines during development.    

**1. Basic figures:**    
This directory contains the scripts that can be used to generate the following figures:     
1.1 Variant types-SMuRF (**Fig. 3a,b** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))    
1.2 ClinVar classification-SMuRF (**Fig. 4a,b** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))    
1.3 Enzymatic domains-SMuRF (**Fig. 5a-d and Supplementary Fig. 14** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))    

**2. Additional analyses:**    
This directory contains the scripts that can be used to generate the following figures:     
2.1 gnomAD-SMuRF (**Fig. 3c,d** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3); see **7.5** bloew)     
2.2 ROC of computational predictors and SMuRF (**Fig. 4e,f** [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))     
2.3 computational predictor score-SMuRF (**Fig. 4g-i** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))     

**3. Patient reports:**    
This directory contains the scripts that can be used to generate the following figures:     
3.1 Disease severity-SMuRF (**Fig. 4c** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))     
3.2 Onset-SMuRF (**Fig. 4d** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))     
3.3 CK values-SMuRF (**Supplementary Fig. 4** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))    

**4. Potential novel discoveries:**    
This directory contains the scripts that can be used to generate the following figures:     
4.1 SynVep scores-SMuRF (**Supplementary Fig. 12** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))     
4.2 Phylop-SMuRF (**Supplementary Fig. 2** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))     

**5. Reclassification:**     
This directory contains the scripts that can be used to classify and/or reclassify variants based on the SMuRF functional scores and confidence scores.    
(**Supplementary Fig. 9,10** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))     

**6. ppVSV enrichment:**     
ppVSV enrichment was evaluated with nanopore sequencing provided by Plasmidusaurus. This directory contains scripts to calculate the functional score based on the nanopore sequencing results.   
(**Fig. 6c,d** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))      

**7. updates:**    
This directory contains the scripts that can be used for following jobs:    
7.1 Adding scores from a new computational tool to the spreadsheet.    
7.2 LoGoFunc-SMuRF (**Supplementary Fig. 13** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))    
7.3 Heatmap for comparing predictors (**Supplementary Fig. 16a,b** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))     
7.4 Enzymatic domains-computational scores & ClinVar-computational scores (**Supplementary Fig. 16c-f** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))       
7.5 gnomAD was updated to v4; we generated the gnomADv4-SMuRF figures. (**Fig. 3c,d** in [BioRxiv](https://www.biorxiv.org/content/10.1101/2023.07.12.548370v3))    

**8. SMuRF score generation:**  
This directory contains the scripts to process the output from DiMSum to generate SMuRF scores:  
8.1 Parsing variant count output from the Gargamel pipeline.  
8.2 Splitting the variant count into respective blocks.  
8.3 Replacing sites with no variant counts with the wild-type value  
8.4 Running DiMSum to combine scores from multiple replicates  
8.5 Formatting the output from DiMSum, perform confidence classification and normalization  

**9. Azrael plots (revised):**  
This directory contains the revised scripts that can be used to generate the following figures:  
9.1 Pair-wise replicate comparison (**Revised Figure S2H**).  
9.2.1 Variant types-SMuRF (**Revised Figures 2A and 2B**).  
9.2.2 ClinVar classification-SMuRF (**Revised Figures 3A, 3B, S3G, and S3H**).  
9.2.3 Enzymatic domains-SMuRF (**Revised Figures 5A-D and S7H**).  
9.3 gnomAD-SMuRF (**Revised Figures 2C, 2D, S3C, and S3D**).  
9.4.1 ROC of computational predictors and SMuRF (**Revised Figures 4A and 4B**).  
9.4.2 computational predictor score-SMuRF (**Revised Figures 4D, 4E, and S4B**).  
9.5 Phylop-SMuRF (**Revised Figures S3E and S3F**).  
9.6.1 Disease severity-SMuRF (**Revised Figures 3C and S3I**).  
9.6.2 Onset-SMuRF (**Revised Figure 3D**).  
9.6.3 CK values-SMuRF (**Revised Figures S3L**).  
9.6.4 ClinVar Disease annotation & the GRASP LGMD consortium figures (**Revised Figures S3J and S3K**).  
9.6.5 Cox proportional hazard test for SMuRF score vs age of onset.  
9.7 SMuRF & computational predictors pair-wise correlation heatmap (**Revised Figure S4A**).  
9.8 A potential SMuRF labeling of variants (**Revised Table 1**).  
9.9 SynVep-SMuRF (**Revised Figure S7F**).  
9.10 LoGoFunc-SMuRF (**Revised Figure S7G**).