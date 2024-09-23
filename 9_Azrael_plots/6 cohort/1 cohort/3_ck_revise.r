#Rscript --vanilla 3_ck_revise.r report_and_SMuRF_high_confidence_dimsum



library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

df = read.csv(args[1],header=T,sep="\t")


df=df[df$CK_IU_L_ave != "na",]
df=df[df$fs1 != "na",]
df=df[df$fs2 != "na",]


r = ggplot(df,aes(x=as.numeric(CK_IU_L_ave), y=log2(2^as.numeric(fs1)+2^as.numeric(fs2)),color=Sex))+ geom_point(size=6)+ geom_smooth(method=lm,se=FALSE, size=1,linetype="dashed") + theme_bw()


r = r + ylab("Log2(2^FS1+2^FS2)") + xlab("CK(IU/L) Ave.") + theme(axis.text.x=element_text(size=12))
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


r = r + theme(legend.position=c(0.85, 0.85),legend.title=element_text(colour="black",size=18))

r = r + guides(fill=guide_legend(title="SEX"))

r = r + theme(legend.background=element_rect(linetype="solid",colour="black"))

# Add a grey rectangle box to show the normal range
r <- r +
  annotate("rect", xmin = 24, xmax = 204, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5)




ggsave(filename=paste(basename(args[1]),".ck_revise.pdf",sep = ""))

dfm=df[df$Sex == "M",]
dff=df[df$Sex == "F",]
cor.test(log2(2^as.numeric(df$fs1)+2^as.numeric(df$fs2)), as.numeric(df$CK_IU_L_ave),method="spearman")
cor.test(log2(2^as.numeric(dfm$fs1)+2^as.numeric(dfm$fs2)), as.numeric(dfm$CK_IU_L_ave),method="spearman")
cor.test(log2(2^as.numeric(dff$fs1)+2^as.numeric(dff$fs2)), as.numeric(dff$CK_IU_L_ave),method="spearman")
