#Rscript --vanilla 2_onset_revise.r report_and_SMuRF_high_confidence_dimsum_NEW


library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

df = read.csv(args[1],header=T,sep="\t")
df=df[df$Onset_m != "na",]
df=df[df$fs1 != "na",]
df=df[df$fs2 != "na",]


r = ggplot(df,aes(x=as.numeric(Onset_m), y=log2(2^as.numeric(fs1)+2^as.numeric(fs2)),color=Sex))+ geom_point(size=6)+ geom_smooth(method=lm,se=FALSE, size=1,linetype="dashed") + theme_bw()


r = r + ylab("Log2(2^FS1+2^FS2)") + xlab("onset (m)") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=1, linetype="dashed", color="black",size=0.5)  + ggtitle("Clivar interpretation - All blocks")

r = r + ggtitle("FKRP well-curated clinical cohorts")


r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black",size=1.5))
r = r + theme(axis.text.y=element_text(size=20,colour="black"),axis.title.y=element_text(size=25,colour="black",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=20,colour="black"),axis.title.x =element_text(size=30,colour="black",vjust=1.5))

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(size=20),plot.subtitle=element_text(size=11))

r=r+ylim(-2.5,2)





r = r + theme(legend.text=element_text(colour="black",size=16))


r = r + theme(legend.position=c(0.9, 0.88),legend.title=element_text(colour="black",size=18))

r = r + guides(fill=guide_legend(title="SEX"))

r = r + theme(legend.background=element_rect(linetype="solid",colour="black"))


ggsave(filename=paste(basename(args[1]),".onset_revise.pdf",sep = ""))


dfm=df[df$Sex == "M",]
dff=df[df$Sex == "F",]
cor.test(log2(2^as.numeric(df$fs1)+2^as.numeric(df$fs2)), as.numeric(df$Onset_m),method="spearman")
cor.test(log2(2^as.numeric(dfm$fs1)+2^as.numeric(dfm$fs2)), as.numeric(dfm$Onset_m),method="spearman")
cor.test(log2(2^as.numeric(dff$fs1)+2^as.numeric(dff$fs2)), as.numeric(dff$Onset_m),method="spearman")
