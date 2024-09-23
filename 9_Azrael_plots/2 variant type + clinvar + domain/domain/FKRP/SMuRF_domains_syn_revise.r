#Rscript --vanilla SMuRF_domains_syn_revise.r FKRP_DiMSum_annotated_NEW_FINAL.tsv

library(ggplot2)
library(ggsignif)

#wilcox.test.default()

args = commandArgs(trailingOnly=TRUE)

all_blocks = read.csv(args[1],header=T,sep="\t")
all_blocks = all_blocks[all_blocks$confidence == "HIGH",]

all_blocks_synonymous = subset(all_blocks, classification == 'Synonymous')
all_blocks_synonymous_stem = subset(all_blocks_synonymous, site<=864)
all_blocks_synonymous_zinc = subset(all_blocks_synonymous, site>864 & site<=954)
all_blocks_synonymous_cat = subset(all_blocks_synonymous, site>954)

all_blocks_synonymous_stem$domain="Stem"
all_blocks_synonymous_zinc$domain="Znf"
all_blocks_synonymous_cat$domain="Catalytic excl. Znf"
all_blocks_synonymous=rbind(all_blocks_synonymous_stem,all_blocks_synonymous_zinc,all_blocks_synonymous_cat)
level<-c("Stem","Znf","Catalytic excl. Znf")

r = ggplot(all_blocks_synonymous, aes(x=factor(domain,level), y=fitness_normalized,color=domain)) + theme_bw()
r = r+scale_color_manual(values=c("mediumpurple","lightyellow3", "deepskyblue2"))
r = r+ theme(legend.position="none")

r = r + geom_violin(lwd=0.6)+geom_signif(comparisons = list(c("Stem", "Znf")),y_position=3.4,color="black")+geom_signif(comparisons = list(c("Stem", "Catalytic excl. Znf")),y_position=3.8,color="black")+geom_signif(comparisons = list(c("Znf", "Catalytic excl. Znf")),y_position=3,color="black")
r = r + geom_boxplot(width=0.1,outlier.size=3)

r = r + scale_y_continuous("SMuRF score") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=0, linetype="dashed", color="black",size=0.5) + scale_x_discrete(name="") + ggtitle("Clivar interpretation - All blocks")

r = r + ggtitle("FKRP Synonymous Variants")


r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black",size=1.5))
r = r + theme(axis.text.y=element_text(size=20,colour="black"),axis.title.y=element_text(size=25,colour="black",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=20,colour="black",hjust=0.6),axis.title.x =element_text(size=10,colour="black"))

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(size=20),plot.subtitle=element_text(size=11))


give.n <- function(x){
  return(c(y = -4.8, label =length(x), size=8))
}


r = r+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))


ggsave(filename="SMuRF_domains_syn_rev_high_confi.pdf")
