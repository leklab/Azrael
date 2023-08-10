#Rscript --vanilla SMuRF_predictors_final_all_test.r gnomad_fkrp-1_annotated_v2.tsv
#Rscript --vanilla SMuRF_predictors_final_all_test.r gnomad_large1_annotated_v2.tsv



library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
scores = read.csv(args[1],header=T,sep="\t")

#1# for EVE
#2# for REVEL
#3# for CADD
#4# for metaSVM
#5# for MVP
#6# for MutScore

#1#
#1#  scores_EVE = subset(scores, EVE!=0)
#1# scores_missense = subset(scores_EVE, classification == 'Missense')

#2346#
scores_missense = subset(scores, classification == 'Missense')

#5#

#5#scores_MVP = subset(scores, mvp_rankscore!=0)
#5#scores_missense = subset(scores_MVP, classification == 'Missense')










#1#
#1# cor.test(scores_missense$functional_score, scores_missense$EVE,method="spearman")

#2#
#2#cor.test(scores_missense$functional_score, scores_missense$revel_score,method="spearman")

#3#
#3# cor.test(scores_missense$functional_score, scores_missense$cadd,method="spearman")

#4#
#4#cor.test(scores_missense$functional_score, scores_missense$meta_svm_score,method="spearman")

#5#
#5#cor.test(scores_missense$functional_score, scores_missense$mvp_rankscore,method="spearman")

#6#
cor.test(scores_missense$functional_score, scores_missense$MutScore,method="spearman")
