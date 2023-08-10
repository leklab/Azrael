#Rscript --vanilla SMuRF_predictors_final.r gnomad_fkrp-1_annotated_v2.tsv
#Rscript --vanilla SMuRF_predictors_final.r gnomad_fkrp-2_annotated_v2.tsv
#Rscript --vanilla SMuRF_predictors_final.r gnomad_large1_annotated_v2.tsv



library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
scores = read.csv(args[1],header=T,sep="\t")

#1# for EVE
#2# for REVEL

#1#
#1# scores_EVE = subset(scores, EVE!=0)
#1# scores_missense = subset(scores_EVE, classification == 'Missense')

#2#
scores_missense = subset(scores, classification == 'Missense')


#1#
#1# s = ggplot(scores_missense, aes(x=EVE, y=log2(functional_score))) + theme_bw()

#2#
s = ggplot(scores_missense, aes(x=revel_score, y=log2(functional_score))) + theme_bw()




#1#
#1# s = s + ylab("log2(Functional score)")+xlab("EVE score")

#2#
s = s + ylab("log2(Functional score)")+xlab("REVEL score")



s = s + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black"))
s = s + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
s = s + theme(axis.text.x=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.x= element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))


s = s + theme(panel.grid.major=element_line(colour="white",size=0.5))
s = s + theme(panel.grid.minor=element_line(colour="white",size=0.5))


variant_count_number=nrow(scores_missense)


########## FKRP or LARGE1 ######
s = s + labs(title="FKRP", subtitle=paste(variant_count_number,"variants"))
#LARGE1# s = s + labs(title="LARGE1", subtitle=paste(variant_count_number,"variants"))



s = s + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))

s = s  + geom_density_2d_filled(contour_var = "count")
s = s + geom_point(alpha=0.05, color="yellow")+ geom_smooth(method=lm,se=FALSE, size=1,linetype="dashed", color="white")

s = s + theme(legend.text=element_text(colour="black",size=10))


s = s + theme(legend.title=element_text(colour="black",size=10))

s = s + guides(fill=guide_legend(title="Count"))

s = s + theme(legend.background=element_rect(linetype="solid",colour="black"))


#1#
#1# ggsave(filename=paste(basename(args[1]),".EVE.png",sep = ""),dpi=300)

#2#
ggsave(filename=paste(basename(args[1]),".REVEL.png",sep = ""),dpi=300)

#1#
#1# cor.test(scores_missense$functional_score, scores_missense$EVE,method="spearman")

#2#
cor.test(scores_missense$functional_score, scores_missense$revel_score,method="spearman")
