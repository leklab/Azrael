#Rscript --vanilla SMuRF_domains_large1_non_revise.r LARGE1_DiMSum_annotated_v2.tsv




library(ggplot2)
library(ggsignif)

#wilcox.test.default()

args = commandArgs(trailingOnly=TRUE)

all_blocks = read.csv(args[1],header=T,sep="\t")
all_blocks = all_blocks[all_blocks$confidence == "HIGH",]




all_blocks_Nonsense = subset(all_blocks, classification == 'Nonsense')
all_blocks_Nonsense_Nterm = subset(all_blocks_Nonsense, site<=411)
all_blocks_Nonsense_Xyl = subset(all_blocks_Nonsense, site>411 & site<=1239)
all_blocks_Nonsense_Glu = subset(all_blocks_Nonsense, site>1239)




all_blocks_Nonsense_Nterm$domain="N-Term"
all_blocks_Nonsense_Xyl$domain="XylT"
all_blocks_Nonsense_Glu$domain="GlcAT"



all_blocks_Nonsense=rbind(all_blocks_Nonsense_Nterm,all_blocks_Nonsense_Xyl,all_blocks_Nonsense_Glu)
level<-c("N-Term","XylT","GlcAT")

###violin
r = ggplot(all_blocks_Nonsense, aes(x=factor(domain,level), y=fitness_normalized,color=domain)) + theme_bw()
r = r + scale_color_manual(values=c("mediumpurple", "lightyellow3","deepskyblue2"))
r = r + theme(legend.position="none")
r = r + geom_violin(lwd=0.6)+geom_signif(comparisons = list(c("N-Term", "XylT")),y_position=4.5,color="black")+geom_signif(comparisons = list(c("N-Term", "GlcAT")),y_position=4.9,color="black")+geom_signif(comparisons = list(c("XylT", "GlcAT")),y_position=4.1,color="black")
r = r + geom_boxplot(width=0.1,outlier.size=3)
r = r + scale_y_continuous("SMuRF score") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=0, linetype="dashed", color="black",size=0.5) + scale_x_discrete(name="") + ggtitle("Clivar interpretation - All blocks")

r = r + ggtitle("LARGE1 Nonsense Variants")


r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black",size=1.5))
r = r + theme(axis.text.y=element_text(size=20,colour="black"),axis.title.y=element_text(size=25,colour="black",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=20,colour="black"),axis.title.x =element_text(size=10,colour="black",vjust=1.5))

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(size=20),plot.subtitle=element_text(size=11))

give.n <- function(x){
  return(c(y = -4.8, label =length(x), size=8))
}


r = r+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))

ggsave(filename="SMuRF_domains_large1_non_revise_high_confi.pdf")
