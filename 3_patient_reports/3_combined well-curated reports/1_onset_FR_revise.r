#Rscript --vanilla 1_onset_FR_revise.r FKRPregstry_report_and_SMuRF_diumsum_NEW report_and_SMuRF_high_confidence_dimsum_NEW




library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

df = read.csv(args[1],header=T,sep="\t")
df2= read.csv(args[2],header=T,sep="\t")

df=subset(df, select=-c(X10MWFDS, NSAD20itemaddition))

df2=subset(df2, select=-c(Sex, Age, Conditions, CK_IU_L_ave, PMID, note))
colnames(df)[4] ="Onset"
df$Onset=as.numeric(df$Onset)*12

colnames(df2)[4] ="Onset"

df=rbind(df,df2)


df=df[df$Onset != "na",]
df=df[df$fs1 != "na",]
df=df[df$fs2 != "na",]
df$Type=ifelse(df$Allele1=="c.826C>A" & df$Allele2=="c.826C>A", "L276I Homozygous", "Other")
df$Type=ifelse(df$Allele1=="c.826C>A" & df$Allele2!="c.826C>A", "L276I Compound Heterozygous", df$Type)
df$Type=ifelse(df$Allele1!="c.826C>A" & df$Allele2=="c.826C>A", "L276I Compound Heterozygous", df$Type)


r = ggplot(df,aes(x=as.numeric(Onset), y=log2(2^as.numeric(fs1)+2^as.numeric(fs2)),color=Type))+ geom_point(size=6) + theme_bw()


r = r + ylab("Log2(2^FS1+2^FS2)") + xlab("onset (month)") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=1, linetype="dashed", color="black",size=0.5)  + ggtitle("Clivar interpretation - All blocks")

r = r + ggtitle("GRASP-LGMD cohort + another 8 well-curated cohorts")


r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black",size=1.5))
r = r + theme(axis.text.y=element_text(size=18,colour="black"),axis.title.y=element_text(size=18,colour="black",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=18,colour="black"),axis.title.x =element_text(size=18,colour="black",vjust=1))

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(size=14),plot.subtitle=element_text(size=11))

r=r+ylim(-2.5,2)





r = r + theme(legend.text=element_text(colour="black",size=10))
r = r + theme(legend.position=c(0.79, 0.3),legend.title=element_text(colour="black",size=10))
r = r + theme(legend.background=element_rect(linetype="solid",colour="black",fill="NA"))


ggsave(filename=paste(basename(args[1]),".onset.pdf",sep = ""),height=3)
