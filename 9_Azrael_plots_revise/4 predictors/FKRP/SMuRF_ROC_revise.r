#Rscript --vanilla SMuRF_ROC_revise.r FKRP_DiMSum_annotated.tsv

### ROC analysis ###
library(pROC)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
df = read.csv(args[1],header=T,sep="\t")
df = df[df$confidence == "HIGH",]

df <- df[df$classification == "Missense",]
df$functional_score <- df$fitness_normalized

clinvar_pos <- c('Pathogenic/Likely pathogenic', 'Pathogenic', 'Likely pathogenic')

# Plot ROC curve

df$response <- ifelse(df$Clinvar_clinical_significance %in% clinvar_pos, 1, 2) # Code cases (Clinvar +) as 1
df$response[df$gnomad_v4_genome_af > 0 & !df$Clinvar_clinical_significance %in% clinvar_pos] <- 0 # Controls as present in gnomAD v4 genome and not Clinvar +
df <- df[!df$response == 2,]

df <- df %>%
  mutate(across(c(meta_svm_score, cadd, revel_score, functional_score, mvp_rankscore, EVE, MutScore, primateAI_3D, Maverick_RecScore, alphamissense, ESM1b_score), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(across(c(meta_svm_score, cadd, revel_score, functional_score, mvp_rankscore, EVE, MutScore, primateAI_3D, Maverick_RecScore, alphamissense, ESM1b_score), as.numeric)) %>%
  filter(meta_svm_score != 0 &
           cadd != 0 &
           revel_score != 0 &
           functional_score != 0 &
           mvp_rankscore != 0 &
           EVE != 0 &
           MutScore != 0 &
           primateAI_3D != 0 &
           Maverick_RecScore != 0 &
           alphamissense != 0 &
           ESM1b_score != 0)

# Compute ROC curves for all predictors
roc_data_metaSVM <- roc(df$response, df$meta_svm_score)
roc_data_CADD <- roc(df$response, df$cadd)
roc_data_REVEL <- roc(df$response, df$revel_score)
roc_data_smurf <- roc(df$response, df$functional_score)
roc_data_gMVP_rankscore <- roc(df$response, df$mvp_rankscore)
roc_data_eve <- roc(df$response, df$EVE)
roc_data_mutscore <- roc(df$response, df$MutScore)
roc_primateai <- roc(df$response,df$primateAI_3D)
roc_maverick <- roc(df$response,as.numeric(df$Maverick_RecScore))
roc_alphamis <- roc(df$response,as.numeric(df$alphamissense))
roc_esm1b <- roc(df$response,df$ESM1b_score)


# Plot ROC curves

png(filename=paste(basename(args[1]),".ROC_revise_high_confidence.png",sep = ""),width=4000, height=4000,res=1200)
par(bg=NA)

plot(roc_data_metaSVM, col = "lightblue",lwd=2,xlab="",ylab="")


plot(roc_data_CADD, col = "gainsboro",lwd=2, add = TRUE)

plot(roc_data_gMVP_rankscore, col = "lavender",lwd=2, add = TRUE)
plot(roc_data_REVEL, col = "green",lwd=2, add = TRUE)
plot(roc_data_eve, col = "blue",lwd=2, add = TRUE)
plot(roc_data_mutscore, col = "lightpink",lwd=2, add = TRUE)

plot(roc_primateai, col = "purple",lwd=2, add = TRUE)
plot(roc_maverick, col = "black",lwd=2, add = TRUE)

plot(roc_alphamis, col = "mistyrose",lwd=2, add = TRUE)

plot(roc_esm1b, col = "darkgreen",lwd=2, add = TRUE)


plot(roc_data_smurf, col = "orange",lwd=2,add = TRUE)





#Add legend with AUC values
legend_text <- c(paste0("MVP (AUC=", round(auc(roc_data_gMVP_rankscore), 2), ")"),
                 paste0("metaSVM (AUC=", round(auc(roc_data_metaSVM), 2), ")"),
                 paste0("ESM1b (AUC=", round(auc(roc_esm1b), 2), ")"),
                 paste0("CADD (AUC=", round(auc(roc_data_CADD), 2), ")"),
                 paste0("MAVERICK (AUC=", round(auc(roc_maverick), 2), ")"),
                 paste0("MutScore (AUC=", round(auc(roc_data_mutscore), 2), ")"),
                 paste0("AlphaMissense (AUC=", round(auc(roc_alphamis), 2), ")"),
                 paste0("PrimateAI_3D (AUC=", round(auc(roc_primateai), 2), ")"),
                 paste0("EVE (AUC=", round(auc(roc_data_eve), 2), ")"),
                 paste0("REVEL (AUC=", round(auc(roc_data_REVEL), 2), ")"),
                 paste0("SMuRF (AUC=", round(auc(roc_data_smurf), 2), ")"))
legend("bottomright", legend = legend_text, col = c("lavender", "lightblue","darkgreen", "gainsboro", "black", 'lightpink','mistyrose',"purple" ,"blue" ,"green",'orange'), lwd = 4,cex=0.34)



##################title##################

title(xlab = "Specificity", line=2.8 , adj=0.5, cex.axis=0.8, font.axis=2, font.lab=2)
title(ylab = "Sensitivity", line=1.9 , adj=0.5 ,cex.axis=0.8, font.axis=2, font.lab=2)

title(main="FKRP", cex.main = 0.8, line=2.2, adj=0.02)
#title(main="LARGE1", cex.main = 0.8, line=2.2, adj=0.02)
