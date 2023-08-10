#Rscript --vanilla 1_bar.r report_and_SMuRF


library(ggplot2)
library(ggsignif)

#wilcox.test.default()

args = commandArgs(trailingOnly=TRUE)

df = read.csv(args[1],header=T,sep="\t")

df=df[df$fs1 != "na",]
df=df[df$fs2 != "na",]



df$Conditions=ifelse(df$Conditions=="MEB", "Severe", df$Conditions)
df$Conditions=ifelse(df$Conditions=="WWS", "Severe", df$Conditions)
df$Conditions=ifelse(df$Conditions=="FCMD", "Severe", df$Conditions)
df$Conditions=ifelse(df$Conditions=="CMD", "Intermediate", df$Conditions)
df$Conditions=ifelse(df$Conditions=="LGMD", "Mild", df$Conditions)

order=c("Severe","Intermediate","Mild")


r = ggplot(df, aes(x=factor(Conditions,order), y=log2(as.numeric(fs1)+as.numeric(fs2)),fill=Sex)) + theme_bw()



r = r + geom_boxplot() + geom_signif(comparisons = list(c("Intermediate", "Mild")),y_position=1.4,tip_length = c(0.05, 0.05))+ geom_signif(comparisons = list(c("Severe", "Mild")),y_position=1,tip_length = c(0.05, 0.05))+ geom_signif(comparisons = list(c("Severe", "Intermediate")),y_position=0.6,tip_length = c(0.05, 0.05))


r = r + ylab("log2(FS1+FS2)") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=1, linetype="dashed", color="blue") + scale_x_discrete(name="") + ggtitle("Clivar interpretation - All blocks")

r = r + ggtitle("FKRP well-curated clinical cohorts")


r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="white"))
r = r + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=14,colour="black",family="Decima Mono Pro"),axis.title.x =element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))
r=r+ylim(-2.5,2)

give.n <- function(x){
  return(c(y = -2.5, label =length(x), size=5))
}


r = r+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))


r = r + theme(legend.text=element_text(colour="black",size=16))


r = r + theme(legend.position=c(0.15, 0.9),legend.title=element_text(colour="black",size=18))

r = r + guides(fill=guide_legend(title="SEX"))

r = r + theme(legend.background=element_rect(linetype="solid",colour="black"))


ggsave(filename=paste(basename(args[1]),".bar.png",sep = ""),dpi=300)
