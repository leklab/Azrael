#Rscript --vanilla 2_onset.r report_and_SMuRF



library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

df = read.csv(args[1],header=T,sep="\t")


df=df[df$Onset_m != "na",]
df=df[df$fs1 != "na",]
df=df[df$fs2 != "na",]


r = ggplot(df,aes(x=as.numeric(Onset_m), y=log2(as.numeric(fs1)+as.numeric(fs2)),color=Sex))+ geom_point(size=6)+ geom_smooth(method=lm,se=FALSE, size=1,linetype="dashed") + theme_bw()


r = r + ylab("log2(FS1+FS2)") + xlab("onset (m)") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=1, linetype="dashed", color="blue")  + ggtitle("Clivar interpretation - All blocks")

r = r + ggtitle("FKRP well-curated clinical cohorts")


r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="white"))
r = r + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=14,colour="black",family="Decima Mono Pro"),axis.title.x =element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))

r=r+ylim(-2.5,2)





r = r + theme(legend.text=element_text(colour="black",size=16))


r = r + theme(legend.position=c(0.85, 0.85),legend.title=element_text(colour="black",size=18))

r = r + guides(fill=guide_legend(title="SEX"))

r = r + theme(legend.background=element_rect(linetype="solid",colour="black"))


ggsave(filename=paste(basename(args[1]),".onset.png",sep = ""),dpi=300)


dfm=df[df$Sex == "M",]
dff=df[df$Sex == "F",]
#cor.test(log2(as.numeric(df$fs1)+as.numeric(df$fs2)), as.numeric(df$Onset_m),method="spearman")
#cor.test(log2(as.numeric(dfm$fs1)+as.numeric(dfm$fs2)), as.numeric(dfm$Onset_m),method="spearman")
#cor.test(log2(as.numeric(dff$fs1)+as.numeric(dff$fs2)), as.numeric(dff$Onset_m),method="spearman")
