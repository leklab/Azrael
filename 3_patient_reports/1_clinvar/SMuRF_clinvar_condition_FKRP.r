#Rscript --vanilla SMuRF_clinvar_condition_FKRP.r FKRP_1_conditions
#Rscript --vanilla SMuRF_clinvar_condition_FKRP.r FKRP_2_conditions


library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

df = read.csv(args[1],header=T,sep="\t")

df$simp=ifelse(df$Clinvar_conditions=="Walker-Warburg congenital muscular dystrophy", "WWS", "other")
df$simp=ifelse(df$Clinvar_conditions=="Autosomal recessive limb-girdle muscular dystrophy" | df$Clinvar_conditions=="Autosomal recessive limb-girdle muscular dystrophy type 2I", "LGMD", df$simp)
df$simp=ifelse(df$Clinvar_conditions=="Muscular dystrophy-dystroglycanopathy type B5" | df$Clinvar_conditions=="Muscular dystrophy-dystroglycanopathy (congenital with brain and eye anomalies), type A5"| df$Clinvar_conditions=="Muscular dystrophy-dystroglycanopathy (congenital with brain and eye anomalies), type A1", "MDDG", df$simp)


order=c("WWS","MDDG","LGMD","other")


r = ggplot(df, aes(x=factor(simp,order), y=log2(functional_score),fill=Clinvar_conditions)) + theme_bw()
#r = r+scale_fill_manual(values=c("mediumpurple","lightyellow", "deepskyblue2"))
r = r+ theme(legend.position="none")
r = r + geom_boxplot()
r = r + scale_y_continuous("log2(Functional score)") + theme(axis.text.x=element_text(size=12))
r = r + geom_hline(yintercept=0, linetype="dashed", color="blue") + scale_x_discrete(name="") + ggtitle("Clivar interpretation - All blocks")

r = r + ggtitle("FKRP ClinVar Conditions (P, P/LP, LP)")


r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="white"))
r = r + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=14,colour="black",family="Decima Mono Pro"),axis.title.x =element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))


give.n <- function(x){
  return(c(y = 0.5, label =length(x), size=5))
}


r = r+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))


ggsave(filename=paste(basename(args[1]),".clinvar_con.png",sep = ""),dpi=300)
