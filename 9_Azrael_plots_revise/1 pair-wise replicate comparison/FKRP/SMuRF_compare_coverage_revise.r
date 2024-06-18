#Rscript --vanilla SMuRF_compare_coverage_revise.r normalized_FKRP-1_func normalized_FKRP-2_func

library(ggplot2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)


exp1 = read.csv(args[1],header=T,sep="\t")
exp2 = read.csv(args[2],header=T,sep="\t")

exp2$b2high_reads=exp2$high_reads


merged <- left_join(exp1, exp2, by = c("site","codon","WT_nt", "Variant_nt", "WT_AA", "Variant_AA"))

exp1$b2high_reads <- merged$b2high_reads


exp1 <- (exp1 %>%
         arrange(factor(classification, levels = c("Missense", "Synonymous", "Nonsense", "Start-loss"))))


s = ggplot(exp1, aes(x=high_reads, y=b2high_reads, color=classification)) + theme_bw()+ scale_color_manual(values = c("Missense" = "lightblue","Nonsense"="red","Start-loss"="darkred","Synonymous"="green"))
s = s + geom_point(size=2)

################
#s = s + labs(caption = paste("Spearman's ρ:", round(cor(exp1$high_reads, exp1$b2high_reads, method = "spearman"), 2)))+theme(plot.caption.position = "plot",plot.caption = element_text(vjust=18,size=17,hjust=0.98))
#s = s + labs(caption = paste("Spearman's ρ:", round(cor(exp1$high_reads, exp1$b2high_reads, method = "spearman",use = "complete.obs"), 2)))+theme(plot.caption.position = "plot",plot.caption = element_text(vjust=18,size=17,hjust=0.16))
s = s + labs(caption = paste("Spearman's ρ:", round(cor(exp1$high_reads, exp1$b2high_reads, method = "spearman",use = "complete.obs"), 2)))+theme(plot.caption.position = "plot",plot.caption = element_text(vjust=18,size=17,hjust=0.95))



#################
#s= s+xlab("FKRP Rep0 High Reads")+ylab("FKRP Rep1 High Reads")
#s= s+xlab("FKRP Rep0 High Reads")+ylab("FKRP Rep2 High Reads")
s= s+xlab("FKRP Rep1 High Reads")+ylab("FKRP Rep2 High Reads")


s = s + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black"))
s = s + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
s = s + theme(axis.text.x=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.x= element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))

######
######s= s + theme(legend.position="none")

s= s + theme(legend.text=element_text(colour="black",family="Decima Mono Pro",size=15))
s= s + theme(legend.position=c(0.19, 0.86),legend.title=element_text(colour="black",family="Decima Mono Pro",size=17))
s= s + theme(legend.background=element_rect(linetype="solid",colour="black",fill="NA"))




s = s + theme(panel.grid.major=element_line(colour="white",size=0.5))
s = s + theme(panel.grid.minor=element_line(colour="white",size=0.5))

#################
#ggsave(filename="compare_coverage_revise_0v1.png",dpi=300)
#ggsave(filename="compare_coverage_revise_0v2.png",dpi=300)
ggsave(filename="compare_coverage_revise_1v2.png",dpi=300)

cor.test(exp1$high_reads, exp1$b2high_reads, method="spearman")
