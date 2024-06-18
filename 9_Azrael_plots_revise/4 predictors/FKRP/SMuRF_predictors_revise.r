#Rscript --vanilla SMuRF_predictors_revise.r FKRP_DiMSum_annotated.tsv


library(ggplot2)


args = commandArgs(trailingOnly=TRUE)
scores = read.csv(args[1],header=T,sep="\t")
scores = scores[scores$confidence == "HIGH",]

#1# for EVE
#2# for REVEL
#3# for primateAI_3D
#4# for MAVERICK
#5# for AlphaMissense

#1#
scores_EVE = subset(scores, EVE!=0)
scores_missense = subset(scores_EVE, classification == 'Missense')

#2345#
#2345#scores_missense = subset(scores, classification == 'Missense')




#1#
s = ggplot(scores_missense, aes(x=EVE, y=fitness_normalized)) + theme_bw()

#2#
#2#s = ggplot(scores_missense, aes(x=revel_score, y=fitness_normalized)) + theme_bw()

#3#
#3#s = ggplot(scores_missense, aes(x=primateAI_3D, y=fitness_normalized)) + theme_bw()

#4#
#4#s = ggplot(scores_missense, aes(x=as.numeric(Maverick_RecScore), y=fitness_normalized)) + theme_bw()

#5#
#5#s = ggplot(scores_missense, aes(x=as.numeric(alphamissense), y=fitness_normalized)) + theme_bw()


#1#
s = s + ylab("SMuRF score")+xlab("EVE score")

#2#
#2#s = s + ylab("SMuRF score")+xlab("REVEL score")

#3#
#3#s = s + ylab("SMuRF score")+xlab("primateAI_3D")

#4#
#4#s = s + ylab("SMuRF score")+xlab("MAVERICK")

#5#
#5#s = s + ylab("SMuRF score")+xlab("AlphaMissense")


s = s + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black"))
s = s + theme(axis.text.y=element_text(size=20,colour="black",family="Atlas Grotesk Web Bold"),axis.title.y=element_text(size=30,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
s = s + theme(axis.text.x=element_text(size=20,colour="black",family="Atlas Grotesk Web Bold"),axis.title.x= element_text(size=30,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))




s = s + theme(panel.grid.major=element_line(colour="white",size=0.5))
s = s + theme(panel.grid.minor=element_line(colour="white",size=0.5))


variant_count_number=nrow(scores_missense)




s = s + labs(title="FKRP", subtitle=paste(variant_count_number,"variants"))


s = s + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=20),plot.subtitle=element_text(family="Atlas Grotesk Web Bold", size=18))

s = s  + geom_density_2d_filled(contour_var = "count")
s = s + geom_point(alpha=0.05, color="yellow")+ geom_smooth(method=lm,se=FALSE, size=1,linetype="dashed", color="white")

s = s + theme(legend.text=element_text(family="Atlas Grotesk Web",colour="black",size=8))+ theme(legend.position="top")


s = s + theme(legend.title=element_text(family="Atlas Grotesk Web",colour="black",size=8))

s = s + guides(fill=guide_legend(title="Count"))

s = s + theme(legend.background=element_rect(linetype="solid",colour="black"))


#1#
ggsave(filename=paste(basename(args[1]),".EVE_high_confi.png",sep = ""),dpi=300, height=10)

#2#
#2#ggsave(filename=paste(basename(args[1]),".REVEL_high_confi.png",sep = ""),dpi=300, height=10)

#3#
#3#ggsave(filename=paste(basename(args[1]),".primateAI_3D_high_confi.png",sep = ""),dpi=300, height=10)

#4#
#4#ggsave(filename=paste(basename(args[1]),".MAVERICK.png",sep = ""),dpi=300, height=10)

#5#
#5#ggsave(filename=paste(basename(args[1]),".AlphaMissense_high_confi.png",sep = ""),dpi=300, height=10)
