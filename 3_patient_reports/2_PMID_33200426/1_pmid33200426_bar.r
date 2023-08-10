#Rscript --vanilla 1_pmid33200426_bar.r PMID_33200426_SMuRF



library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

df = read.csv(args[1],header=T,sep="\t")



order=c("CMD-MR","CMD","LGMD-MR","LGMD")


r = ggplot(df, aes(x=factor(Conditions,order), y=log2(fs1+fs2),fill=Sex)) + theme_bw()
#r = r+scale_fill_manual(values=c("mediumpurple","lightyellow", "deepskyblue2"))
r = r+ theme(legend.position="none")
r = r + geom_boxplot()
r = r + scale_y_continuous("log2(FS1+FS2)") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=1, linetype="dashed", color="blue") + scale_x_discrete(name="") + ggtitle("Clivar interpretation - All blocks")

r = r + ggtitle("FKRP PMID 33200426")


r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="white"))
r = r + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=14,colour="black",family="Decima Mono Pro"),axis.title.x =element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))
r=r+ylim(-2,2)

give.n <- function(x){
  return(c(y = -2, label =length(x), size=5))
}


r = r+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))


r = r + theme(legend.text=element_text(colour="black",size=16))


r = r + theme(legend.position=c(0.15, 0.85),legend.title=element_text(colour="black",size=18))

r = r + guides(fill=guide_legend(title="SEX"))

r = r + theme(legend.background=element_rect(linetype="solid",colour="black"))


ggsave(filename=paste(basename(args[1]),".33200426.png",sep = ""),dpi=300)
