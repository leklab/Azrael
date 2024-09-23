#Rscript --vanilla SMuRF_LARGE1_combine_plot_revise.r LARGE1_DiMSum_annotated.tsv

library(ggplot2)
library(ggsignif)

#wilcox.test.default()

args = commandArgs(trailingOnly=TRUE)

all_blocks = read.csv(args[1],header=T,sep="\t")
#all_blocks = all_blocks[all_blocks$confidence == "HIGH",]

level_order <- c('Synonymous', 'Missense', 'Nonsense', 'Start-loss')


p = ggplot(all_blocks, aes(x=factor(classification,level_order), y=fitness_normalized,color=classification)) + theme_bw()

p = p + geom_violin(lwd=0.6) + geom_signif(comparisons = list(c("Synonymous", "Missense")),y_position=3,color="black")+geom_signif(comparisons = list(c("Synonymous", "Nonsense")),y_position=3.4,color="black")+geom_signif(comparisons = list(c("Synonymous", "Start-loss")),y_position=3.8,color="black")+geom_signif(comparisons = list(c("Nonsense", "Start-loss")),y_position=3,color="black")
p=p+geom_boxplot(width=0.1,outlier.size=3)


p = p+scale_color_manual(values=c("lightblue", "red","darkred","darkgreen"))
p = p+ theme(legend.position="none")


p = p + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black",size=1.5))
p = p + theme(axis.text.y=element_text(size=20,colour="black"),axis.title.y=element_text(size=25,colour="black",vjust=1.5))
p = p + theme(axis.text.x=element_text(size=17,colour="black"),axis.title.x = element_blank())

p = p + theme(panel.grid.major=element_line(colour="white",size=0.5))
p = p + theme(panel.grid.minor=element_line(colour="white",size=0.5))
p = p + theme(axis.ticks=element_blank(),plot.title=element_text(size=14),plot.subtitle=element_text(size=11))



variant_count_number=nrow(all_blocks)
p = p + labs(title="LARGE1 - All Blocks", subtitle=paste(variant_count_number,"variants"))
p = p + ylab('SMuRF score') + geom_hline(yintercept=0, linetype="dashed",color="black",size=0.5)
p= p + theme(plot.title=element_text(size=20))
p= p + theme(plot.subtitle=element_text(size=18))

give.n <- function(x){
  return(c(y = -4.8, label =length(x),size=8))
}

p=p+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))


#ggsave(filename="SMuRF_large1_scores_combined_style_withnumber_revised_high_confidence.pdf")
ggsave(filename="SMuRF_large1_scores_combined_style_withnumber_revised.pdf")


###
Synonymous_score=all_blocks[all_blocks$classification=="Synonymous",]
Missense_score=all_blocks[all_blocks$classification=="Missense",]
Nonsense_score=all_blocks[all_blocks$classification=="Nonsense",]
Startloss_score=all_blocks[all_blocks$classification=="Start-loss",]


print(c("synmedian:", median(Synonymous_score$fitness_normalized)))
synmedian = apply(matrix(sample(Synonymous_score$fitness_normalized, rep=TRUE, 10^4*length(Synonymous_score$fitness_normalized)), nrow=10^4), 1, median)
quantile(synmedian, c(.025, 0.975))

print(c("nonmedian:", median(Nonsense_score$fitness_normalized)))
nonmedian = apply(matrix(sample(Nonsense_score$fitness_normalized, rep=TRUE, 10^4*length(Nonsense_score$fitness_normalized)), nrow=10^4), 1, median)
quantile(nonmedian, c(.025, 0.975))

print(c("stamedian:", median(Startloss_score$fitness_normalized)))
stamedian = apply(matrix(sample(Startloss_score$fitness_normalized, rep=TRUE, 10^4*length(Startloss_score$fitness_normalized)), nrow=10^4), 1, median)
quantile(stamedian, c(.025, 0.975))
