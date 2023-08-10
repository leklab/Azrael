#Rscript --vanilla Compare_LARGE1.r normalized_LARGE1_block1_func normalized_LARGE1_block2_func normalized_LARGE1_block3_func normalized_LARGE1_block4_func normalized_LARGE1_block5_func normalized_LARGE1_block6_func normalized_LARGE1_block7_func normalized_LARGE1_block8_func normalized_LARGE1_block9_func normalized_LARGE1_block10_func old_normalized_LARGE1_block1_func old_normalized_LARGE1_block2_func old_normalized_LARGE1_block3_func old_normalized_LARGE1_block4_func old_normalized_LARGE1_block5_func old_normalized_LARGE1_block6_func old_normalized_LARGE1_block7_func old_normalized_LARGE1_block8_func old_normalized_LARGE1_block9_func old_normalized_LARGE1_block10_func

library(ggplot2)
library(dplyr)

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

old_block1 = read.csv(args[11],header=T,sep="\t")
old_block2 = read.csv(args[12],header=T,sep="\t")
old_block3 = read.csv(args[13],header=T,sep="\t")
old_block4 = read.csv(args[14],header=T,sep="\t")
old_block5 = read.csv(args[15],header=T,sep="\t")
old_block6 = read.csv(args[16],header=T,sep="\t")
old_block7 = read.csv(args[17],header=T,sep="\t")
old_block8 = read.csv(args[18],header=T,sep="\t")
old_block9 = read.csv(args[19],header=T,sep="\t")
old_block10 = read.csv(args[20],header=T,sep="\t")



new=rbind(block1,block2,block3,block4,block5,block6,block7,block8,block9,block10)
old=rbind(old_block1,old_block2,old_block3,old_block4,old_block5,old_block6,old_block7,old_block8,old_block9,old_block10)
old$old_functional_score=old$functional_score



merged <- left_join(new, old, by = c("site","codon","WT_nt", "Variant_nt", "WT_AA", "Variant_AA"))

new$old_functional_score <- merged$old_functional_score




s = ggplot(new, aes(x=log2(functional_score), y=log2(old_functional_score))) + theme_bw()
s = s + geom_point(alpha=0.65)
s = s + xlab("New log2(Functional score)")+ylab("Old log2(Functional score)")



s = s + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black"))
s = s + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
s = s + theme(axis.text.x=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.x= element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))


s = s + theme(panel.grid.major=element_line(colour="white",size=0.5))
s = s + theme(panel.grid.minor=element_line(colour="white",size=0.5))
ggsave(filename="compare_large1.png",dpi=300)

cor.test(new$functional_score, new$old_functional_score, method="spearman")
