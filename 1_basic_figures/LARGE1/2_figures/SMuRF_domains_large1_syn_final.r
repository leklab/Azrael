#Rscript --vanilla SMuRF_domains_large1_syn_final.r normalized_LARGE1_block1_func normalized_LARGE1_block2_func normalized_LARGE1_block3_func normalized_LARGE1_block4_func normalized_LARGE1_block5_func normalized_LARGE1_block6_func normalized_LARGE1_block7_func normalized_LARGE1_block8_func normalized_LARGE1_block9_func normalized_LARGE1_block10_func


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
block7 = read.csv(args[7],header=T,sep="\t")
block8 = read.csv(args[8],header=T,sep="\t")
block9 = read.csv(args[9],header=T,sep="\t")
block10 = read.csv(args[10],header=T,sep="\t")


block1$block = 'block1'
block2$block = 'block2'
block3$block = 'block3'
block4$block = 'block4'
block5$block = 'block5'
block6$block = 'block6'
block7$block = 'block7'
block8$block = 'block8'
block9$block = 'block9'
block10$block = 'block10'



all_blocks = rbind(block1,block2,block3,block4,block5,block6,block7,block8,block9,block10)




all_blocks_synonymous = subset(all_blocks, classification == 'Synonymous')
all_blocks_synonymous_Nterm = subset(all_blocks_synonymous, site<=411)
all_blocks_synonymous_Xyl = subset(all_blocks_synonymous, site>411 & site<=1239)
all_blocks_synonymous_Glu = subset(all_blocks_synonymous, site>1239)




all_blocks_synonymous_Nterm$domain="N-Term"
all_blocks_synonymous_Xyl$domain="XylT"
all_blocks_synonymous_Glu$domain="GlcAT"



all_blocks_synonymous=rbind(all_blocks_synonymous_Nterm,all_blocks_synonymous_Xyl,all_blocks_synonymous_Glu)
level<-c("N-Term","XylT","GlcAT")

r = ggplot(all_blocks_synonymous, aes(x=factor(domain,level), y=log2(functional_score),fill=domain)) + theme_bw()
r = r+scale_fill_manual(values=c("mediumpurple", "lightyellow","deepskyblue2"))
r = r+ theme(legend.position="none")
r = r + geom_boxplot()+geom_signif(comparisons = list(c("N-Term", "XylT")),y_position=4.5)+geom_signif(comparisons = list(c("N-Term", "GlcAT")),y_position=4.9)+geom_signif(comparisons = list(c("XylT", "GlcAT")),y_position=4.1)
r = r + scale_y_continuous("log2(Functional score)") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=0, linetype="dashed", color="blue") + scale_x_discrete(name="") + ggtitle("Clivar interpretation - All blocks")

r = r + ggtitle("LARGE1 Synonymous Variants")


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

ggsave(filename="SMuRF_domains_large1_syn.png",dpi=300)
