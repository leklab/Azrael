#Rscript --vanilla SMuRF_domains_mis_final.r normalized_FKRP_block1_func normalized_FKRP_block2_func  normalized_FKRP_block3_func normalized_FKRP_block4_func normalized_FKRP_block5_func normalized_FKRP_block6_func


library(ggplot2)
library(ggsignif)

#wilcox.test.default()

args = commandArgs(trailingOnly=TRUE)

block1 = read.csv(args[1],header=T,sep="\t")
block2 = read.csv(args[2],header=T,sep="\t")
block3 = read.csv(args[3],header=T,sep="\t")
block4 = read.csv(args[4],header=T,sep="\t")
block5 = read.csv(args[5],header=T,sep="\t")
block6 = read.csv(args[6],header=T,sep="\t")


block1$block = 'block1'
block2$block = 'block2'
block3$block = 'block3'
block4$block = 'block4'
block5$block = 'block5'
block6$block = 'block6'


all_blocks = rbind(block1,block2,block3,block4,block5,block6)




all_blocks_missense = subset(all_blocks, classification == 'Missense')
all_blocks_missense_stem = subset(all_blocks_missense, site<=864)
all_blocks_missense_zinc = subset(all_blocks_missense, site>864 & site<=954)
all_blocks_missense_cat = subset(all_blocks_missense, site>954)


all_blocks_missense_stem$domain="Stem"
all_blocks_missense_zinc$domain="Znf"
all_blocks_missense_cat$domain="Catalytic excl. Znf"
all_blocks_missense=rbind(all_blocks_missense_stem,all_blocks_missense_zinc,all_blocks_missense_cat)
level<-c("Stem","Znf","Catalytic excl. Znf")

r = ggplot(all_blocks_missense, aes(x=factor(domain,level), y=log2(functional_score),fill=domain)) + theme_bw()
r = r+scale_fill_manual(values=c("mediumpurple","lightyellow", "deepskyblue2"))
r = r+ theme(legend.position="none")

r = r + geom_boxplot()+geom_signif(comparisons = list(c("Stem", "Znf")),y_position=4)+geom_signif(comparisons = list(c("Stem", "Catalytic excl. Znf")),y_position=4.4)+geom_signif(comparisons = list(c("Znf", "Catalytic excl. Znf")),y_position=3.6)


r = r + scale_y_continuous("log2(Functional score)") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=0, linetype="dashed", color="blue") + scale_x_discrete(name="") + ggtitle("Clivar interpretation - All blocks")

r = r + ggtitle("FKRP Missense Variants")


r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="white"))
r = r + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=14,colour="black",family="Decima Mono Pro"),axis.title.x =element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))


give.n <- function(x){
  return(c(y = -4.8, label =length(x), size=5))
}


r = r+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))


ggsave(filename="SMuRF_domains_mis.png",dpi=300)
