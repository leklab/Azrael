#Rscript --vanilla SMuRF_LARGE1_combine_plot_final.r normalized_LARGE1_block1_func normalized_LARGE1_block2_func normalized_LARGE1_block3_func normalized_LARGE1_block4_func normalized_LARGE1_block5_func normalized_LARGE1_block6_func normalized_LARGE1_block7_func normalized_LARGE1_block8_func normalized_LARGE1_block9_func normalized_LARGE1_block10_func


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




level_order <- c('Synonymous', 'Missense', 'Nonsense', 'Start-loss')


p = ggplot(all_blocks, aes(x=factor(classification,level_order), y=log2(functional_score),fill=classification)) + theme_bw()
p = p + geom_boxplot() + geom_signif(comparisons = list(c("Synonymous", "Missense")),y_position=4.1)+geom_signif(comparisons = list(c("Synonymous", "Nonsense")),y_position=4.5)+geom_signif(comparisons = list(c("Synonymous", "Start-loss")),y_position=4.9)+geom_signif(comparisons = list(c("Nonsense", "Start-loss")),y_position=4.1)
p = p+scale_fill_manual(values=c("lightblue", "red","darkred","darkgreen"))
p = p+ theme(legend.position="none")


p = p + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="white"))
p = p + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
p = p + theme(axis.text.x=element_text(size=16,colour="black",family="Decima Mono Pro"),axis.title.x = element_blank())

p = p + theme(panel.grid.major=element_line(colour="white",size=0.5))
p = p + theme(panel.grid.minor=element_line(colour="white",size=0.5))
p = p + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))




variant_count_number=length(readLines(args[1]))+length(readLines(args[2]))+length(readLines(args[3]))+length(readLines(args[4]))+length(readLines(args[5]))+length(readLines(args[6]))+length(readLines(args[7]))+length(readLines(args[8]))+length(readLines(args[9]))+length(readLines(args[10]))-10
p = p + labs(title="LARGE1 - All Blocks", subtitle=paste(variant_count_number,"variants"))
p = p + ylab('log2(Functional score)') + geom_hline(yintercept=0, linetype="dashed",color="black")
p= p + theme(plot.title=element_text(colour="black"))
p= p + theme(plot.subtitle=element_text(colour="black"))

give.n <- function(x){
  return(c(y = -4.8, label =length(x)))
}

p=p+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))


ggsave(filename="SMuRF_large1_scores_combined_style_withnumber.png",dpi=300)

Synonymous_score=all_blocks[all_blocks$classification=="Synonymous",]
Missense_score=all_blocks[all_blocks$classification=="Missense",]
Nonsense_score=all_blocks[all_blocks$classification=="Nonsense",]
Startloss_score=all_blocks[all_blocks$classification=="Start-loss",]


print(c("synmedian:", median(log2(Synonymous_score$functional_score))))
synmedian = apply(matrix(sample(log2(Synonymous_score$functional_score), rep=TRUE, 10^4*length(Synonymous_score$functional_score)), nrow=10^4), 1, median)
quantile(synmedian, c(.025, 0.975))

print(c("nonmedian:", median(log2(Nonsense_score$functional_score))))
nonmedian = apply(matrix(sample(log2(Nonsense_score$functional_score), rep=TRUE, 10^4*length(Nonsense_score$functional_score)), nrow=10^4), 1, median)
quantile(nonmedian, c(.025, 0.975))

print(c("stamedian:", median(log2(Startloss_score$functional_score))))
stamedian = apply(matrix(sample(log2(Startloss_score$functional_score), rep=TRUE, 10^4*length(Startloss_score$functional_score)), nrow=10^4), 1, median)
quantile(stamedian, c(.025, 0.975))
