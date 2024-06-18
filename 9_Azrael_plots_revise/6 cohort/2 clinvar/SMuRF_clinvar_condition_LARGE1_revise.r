#Rscript --vanilla SMuRF_clinvar_condition_LARGE1_revise.r LARGE1_conditions



library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

df = read.csv(args[1],header=T,sep="\t")

df$simp=ifelse(df$Clinvar_conditions=="Retinitis pigmentosa", "RP", "other")
df$simp=ifelse(df$Clinvar_conditions=="Muscular dystrophy-dystroglycanopathy type B6" , "MDDGB6", df$simp)
df$simp=ifelse(df$Clinvar_conditions=="Muscular dystrophy-dystroglycanopathy (congenital with brain and eye anomalies), type A6" , "MDDGA6", df$simp)


order=c("MDDGA6","MDDGB6","RP")


r = ggplot(df, aes(x=factor(simp,order), y=fitness_normalized,fill=Clinvar_conditions)) + theme_bw()
#r = r+scale_fill_manual(values=c("mediumpurple","lightyellow", "deepskyblue2"))
r = r+ theme(legend.position="none")
r = r + geom_boxplot()
r = r + scale_y_continuous("SMuRF score") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=0, linetype="dashed", color="black",size=0.5) + scale_x_discrete(name="") + ggtitle("Clivar interpretation - All blocks")

r = r + ggtitle("LARGE1 ClinVar Conditions (P, P/LP, LP)")


r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black",size=1.5))
r = r + theme(axis.text.y=element_text(size=20,colour="black",family="Atlas Grotesk Web Bold"),axis.title.y=element_text(size=25,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=25,colour="black",family="Atlas Grotesk Web Bold"),axis.title.x =element_text(size=10,colour="black",family="Atlas Grotesk Web Bold",vjust=1))

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))


give.n <- function(x){
  return(c(y = 0.5, label =length(x), size=8))
}


r = r+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))


ggsave(filename=paste(basename(args[1]),".clinvar_con_revise.png",sep = ""),dpi=300,width=6.5,height=14)
