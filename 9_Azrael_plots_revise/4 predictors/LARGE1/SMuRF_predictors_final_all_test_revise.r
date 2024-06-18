#Rscript --vanilla SMuRF_predictors_final_all_test_revise.r LARGE1_DiMSum_annotated_v6.tsv

library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
scores = read.csv(args[1],header=T,sep="\t")
scores = scores[scores$confidence == "HIGH",]

#1# for EVE
#2# for REVEL
#3# for CADD
#4# for metaSVM
#5# for MVP
#6# for MutScore
#7# for PrimateAI_3D
#8# for MAVERICK
#9# for AlphaMissense
#10# for ESM1b

#1#
#1#scores_EVE = subset(scores, EVE!=0)
#1#scores_missense = subset(scores_EVE, classification == 'Missense')

#234678910#
scores_missense = subset(scores, classification == 'Missense')

#5#

#5#scores_MVP = subset(scores, mvp_rankscore!=0)
#5#scores_missense = subset(scores_MVP, classification == 'Missense')




#1#
#1#cor.test(scores_missense$fitness_normalized, scores_missense$EVE,method="spearman")

#2#
cor.test(scores_missense$fitness_normalized, scores_missense$revel_score,method="spearman")

#3#
cor.test(scores_missense$fitness_normalized, scores_missense$cadd,method="spearman")

#4#
cor.test(scores_missense$fitness_normalized, scores_missense$meta_svm_score,method="spearman")

#5#
#5#cor.test(scores_missense$fitness_normalized, scores_missense$mvp_rankscore,method="spearman")

#6#
cor.test(scores_missense$fitness_normalized, scores_missense$MutScore,method="spearman")

#7#
cor.test(scores_missense$fitness_normalized, scores_missense$primateAI_3D,method="spearman")

#8#
cor.test(scores_missense$fitness_normalized, as.numeric(scores_missense$Maverick_RecScore),method="spearman")

#9#
cor.test(scores_missense$fitness_normalized, as.numeric(scores_missense$alphamissense),method="spearman")

#10#
cor.test(scores_missense$fitness_normalized, scores_missense$ESM1b_score ,method="spearman")
