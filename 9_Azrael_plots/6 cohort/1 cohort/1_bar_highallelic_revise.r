#Rscript --vanilla 1_bar_highallelic_revise.r report_and_SMuRF_high_confidence_dimsum

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
df$FS=ifelse(as.numeric(df$fs1)>=as.numeric(df$fs2), as.numeric(df$fs1), as.numeric(df$fs2))
order=c("Severe","Intermediate","Mild")


r = ggplot(df, aes(x=factor(Conditions,order),y=FS,color=Sex)) + theme_bw()



r = r + geom_violin(lwd=1) + geom_signif(comparisons = list(c("Intermediate", "Mild")),y_position=0.6,tip_length = c(0.05, 0.05),size=1,textsize=7,color="black")+ geom_signif(comparisons = list(c("Severe", "Mild")),y_position=0.2,tip_length = c(0.05, 0.05),size=1,textsize=7,color="black")+ geom_signif(comparisons = list(c("Severe", "Intermediate")),y_position=-0.4,tip_length = c(0.05, 0.05),size=1,textsize=7,color="black")
r = r + geom_dotplot(fill="white",binaxis='y', stackdir='center',position=position_dodge(0.9),dotsize=0.5,stackratio=0.4)


r = r + ylab("Max(FS1,FS2)") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=0, linetype="dashed", color="black", size=0.5) + scale_x_discrete(name="")

r = r + ggtitle("FKRP well-curated clinical cohorts")


r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black",size=1.5))
r = r + theme(axis.text.y=element_text(size=20,colour="black"),axis.title.y=element_text(size=25,colour="black",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=25,colour="black"),axis.title.x =element_text(size=10,colour="black",vjust=1))

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(size=20),plot.subtitle=element_text(size=11))
r=r+ylim(-3.5,1)

give.n <- function(x){
  return(c(y = -3.5, label =length(x), size=10))
}


r = r+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))


r = r + theme(legend.text=element_text(colour="black",size=16))


r = r + theme(legend.position=c(0.1, 0.88),legend.title=element_text(colour="black",size=18))

r = r + guides(color=guide_legend(title="SEX",override.aes = aes(label = "")))

r = r + theme(legend.background=element_rect(linetype="solid",colour="black"))


ggsave(filename=paste(basename(args[1]),".bar.revise_high_allele.pdf",sep = ""))
