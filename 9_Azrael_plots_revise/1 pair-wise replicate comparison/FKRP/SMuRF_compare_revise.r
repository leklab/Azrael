#Rscript --vanilla SMuRF_compare_revise.r normalized_FKRP-0_func normalized_FKRP-1_func
#Rscript --vanilla SMuRF_compare_revise.r normalized_FKRP-0_func normalized_FKRP-2_func
#Rscript --vanilla SMuRF_compare_revise.r normalized_FKRP-1_func normalized_FKRP-2_func

library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)


exp1 = read.csv(args[1],header=T,sep="\t")
exp2 = read.csv(args[2],header=T,sep="\t")

exp2$b2functional_score=exp2$functional_score


merged <- left_join(exp1, exp2, by = c("site","codon","WT_nt", "Variant_nt", "WT_AA", "Variant_AA"))

exp1$b2functional_score <- merged$b2functional_score


exp1 <- (exp1 %>%
         arrange(factor(classification, levels = c("Missense", "Synonymous", "Nonsense", "Start-loss"))))


s = ggplot(exp1, aes(x=log2(functional_score), y=log2(b2functional_score), color=classification)) + theme_bw()+ scale_color_manual(values = c("Missense" = "lightblue","Nonsense"="red","Start-loss"="darkred","Synonymous"="green"))
s = s + geom_point(size=2)

################
#s = s + labs(caption = paste("Spearman's ρ:", round(cor(exp1$functional_score, exp1$b2functional_score, method = "spearman"), 2)))+theme(plot.caption.position = "plot",plot.caption = element_text(vjust=18,size=17,hjust=0.98))
#s = s + labs(caption = paste("Spearman's ρ:", round(cor(exp1$functional_score, exp1$b2functional_score, method = "spearman",use = "complete.obs"), 2)))+theme(plot.caption.position = "plot",plot.caption = element_text(vjust=18,size=17,hjust=0.16))
s = s + labs(caption = paste("Spearman's ρ:", round(cor(exp1$functional_score, exp1$b2functional_score, method = "spearman",use = "complete.obs"), 2)))+theme(plot.caption.position = "plot",plot.caption = element_text(vjust=18,size=17,hjust=0.16))



#################
#s= s+xlab("FKRP Rep0 log2(Functional score)")+ylab("FKRP Rep1 log2(Functional score)")
#s= s+xlab("FKRP Rep0 log2(Functional score)")+ylab("FKRP Rep2 log2(Functional score)")
s= s+xlab("FKRP Rep1 log2(Functional score)")+ylab("FKRP Rep2 log2(Functional score)")


s = s + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black"))
s = s + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
s = s + theme(axis.text.x=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.x= element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))

######
s= s + theme(legend.position="none")

#s= s + theme(legend.text=element_text(colour="black",family="Decima Mono Pro",size=15))
#s= s + theme(legend.position=c(0.19, 0.86),legend.title=element_text(colour="black",family="Decima Mono Pro",size=17))
#s= s + theme(legend.background=element_rect(linetype="solid",colour="black",fill="NA"))




s = s + theme(panel.grid.major=element_line(colour="white",size=0.5))
s = s + theme(panel.grid.minor=element_line(colour="white",size=0.5))

#################
#ggsave(filename="compare_revise_0v1.png",dpi=300)
#ggsave(filename="compare_revise_0v2.png",dpi=300)
ggsave(filename="compare_revise_1v2.png",dpi=300)

cor.test(exp1$functional_score, exp1$b2functional_score, method="spearman")
