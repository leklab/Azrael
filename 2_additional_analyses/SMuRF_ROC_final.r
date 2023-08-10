#Rscript --vanilla SMuRF_ROC_final.r gnomad_fkrp-1_annotated_v2.tsv
#Rscript --vanilla SMuRF_ROC_final.r gnomad_fkrp-2_annotated_v2.tsv
#Rscript --vanilla SMuRF_ROC_final.r gnomad_large1_annotated_v2.tsv

### ROC analysis ###
library(pROC)


args = commandArgs(trailingOnly=TRUE)
df = read.csv(args[1],header=T,sep="\t")
df <- df[df$classification == "Missense",]

clinvar_pos <- c('Pathogenic/Likely pathogenic', 'Pathogenic', 'Likely pathogenic')

# Plot ROC curve

df$response <- ifelse(df$Clinvar_clinical_significance %in% clinvar_pos, 1, 0)   #Code cases (Clinvar +) as 1, and controls as 0

# Compute ROC curves for all predictors
roc_data_metaSVM <- roc(df$response, df$meta_svm_score)
roc_data_CADD <- roc(df$response, df$cadd)
roc_data_REVEL <- roc(df$response, df$revel_score)
roc_data_smurf <- roc(df$response, df$functional_score)


dfMVP <- df[df$mvp_rankscore != 0,]
dfEVE <- df[df$EVE != 0,]

roc_data_gMVP_rankscore <- roc(dfMVP$response, dfMVP$mvp_rankscore)
roc_data_eve <- roc(dfEVE$response, dfEVE$EVE)

roc_data_mutscore <- roc(df$response, df$MutScore)


# Plot ROC curves

png(filename=paste(basename(args[1]),".ROC.png",sep = ""),width=4000, height=4000,res=600)

plot(roc_data_metaSVM, col = "lightblue",lwd=4, xlab = "Specificity", ylab = "Sensitivity",cex.lab=2, cex.axis=1.5)
plot(roc_data_CADD, col = "gainsboro",lwd=4, add = TRUE)

plot(roc_data_gMVP_rankscore, col = "lavender",lwd=4, add = TRUE)
plot(roc_data_REVEL, col = "green",lwd=4, add = TRUE)
plot(roc_data_eve, col = "blue",lwd=4, add = TRUE)
plot(roc_data_mutscore, col = "lightpink",lwd=4, add = TRUE)

plot(roc_data_smurf, col = "orange",lwd=4,add = TRUE)



# Add legend with AUC values
legend_text <- c(paste0("metaSVM (AUC=", round(auc(roc_data_metaSVM), 2), ")"),
                 paste0("CADD (AUC=", round(auc(roc_data_CADD), 2), ")"),
                 paste0("MVP (AUC=", round(auc(roc_data_gMVP_rankscore), 2), ")"),
                 paste0("REVEL (AUC=", round(auc(roc_data_REVEL), 2), ")"),
                 paste0("EVE (AUC=", round(auc(roc_data_eve), 2), ")"),
                 paste0("MutScore (AUC=", round(auc(roc_data_mutscore), 2), ")"),
                 paste0("SmuRF (AUC=", round(auc(roc_data_smurf), 2), ")"))
legend("bottomright", legend = legend_text, col = c("lightblue", "gainsboro", "lavender", "green", 'blue','lightpink', 'orange'), lwd = 4,cex=1.25)

#title(main="FKRP", cex.main = 1.2, line=0.6, adj=0.02)
title(main="LARGE1", cex.main = 1.2, line=0.6, adj=0.02)
