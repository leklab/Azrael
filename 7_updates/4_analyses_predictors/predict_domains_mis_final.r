#Rscript --vanilla predict_domains_mis_final.r Supplementary\ Table\ 4_updated_v4.tsv
#Rscript --vanilla predict_domains_mis_final.r Supplementary\ Table\ 5_updated_v5.tsv


library(ggplot2)
library(ggsignif)

#wilcox.test.default()

args = commandArgs(trailingOnly=TRUE)

all_blocks  = read.csv(args[1],header=T,sep="\t")
all_blocks  <- subset(all_blocks , am_pathogenicity != "na")
all_blocks_missense = subset(all_blocks, classification == 'Missense')



# FKRP-specific

# all_blocks_missense_stem = subset(all_blocks_missense, site<=864)
# all_blocks_missense_zinc = subset(all_blocks_missense, site>864 & site<=954)
# all_blocks_missense_cat = subset(all_blocks_missense, site>954)
# all_blocks_missense_stem$domain="Stem"
# all_blocks_missense_zinc$domain="Znf"
# all_blocks_missense_cat$domain="Catalytic excl. Znf"
# all_blocks_missense=rbind(all_blocks_missense_stem,all_blocks_missense_zinc,all_blocks_missense_cat)
# level<-c("Stem","Znf","Catalytic excl. Znf")

# FKRP-specific
# LARGE1-specific

all_blocks_missense_Nterm = subset(all_blocks_missense, site<=411)
all_blocks_missense_Xyl = subset(all_blocks_missense, site>411 & site<=1239)
all_blocks_missense_Glu = subset(all_blocks_missense, site>1239)
all_blocks_missense_Nterm$domain="N-Term"
all_blocks_missense_Xyl$domain="XylT"
all_blocks_missense_Glu$domain="GlcAT"
all_blocks_missense=rbind(all_blocks_missense_Nterm,all_blocks_missense_Xyl,all_blocks_missense_Glu)
level<-c("N-Term","XylT","GlcAT")

# LARGE1-specific


r = ggplot(all_blocks_missense, aes(x=factor(domain,level), y=as.numeric(am_pathogenicity),fill=domain)) + theme_bw()
r = r+scale_fill_manual(values=c("mediumpurple","lightyellow", "deepskyblue2"))
r = r+ theme(legend.position="none")

# FKRP-specific

# r = r + geom_boxplot(lwd=1)+geom_signif(comparisons = list(c("Stem", "Znf")),y_position=1.1,size=1,textsize=5)+geom_signif(comparisons = list(c("Stem", "Catalytic excl. Znf")),y_position=1.2,size=1,textsize=5)+geom_signif(comparisons = list(c("Znf", "Catalytic excl. Znf")),y_position=1.3,size=1,textsize=5)
# r = r + ggtitle("FKRP Missense Variants")

# FKRP-specific

# LARGE1-specific

r = r + geom_boxplot(lwd=1)+geom_signif(comparisons = list(c("N-Term", "XylT")),y_position=1.1,size=1,textsize=5)+geom_signif(comparisons = list(c("N-Term", "GlcAT")),y_position=1.2,size=1,textsize=5)+geom_signif(comparisons = list(c("XylT", "GlcAT")),y_position=1.3,size=1,textsize=5)
r = r + ggtitle("LARGE1 Missense Variants")

# LARGE1-specific

r = r + scale_y_continuous("AlphaMissense") + theme(axis.text.x=element_text(size=12))




r = r + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black",size=1.5))
r = r + theme(axis.text.y=element_text(size=20,colour="black",family="Atlas Grotesk Web Bold"),axis.title.y=element_text(size=25,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
r = r + theme(axis.text.x=element_text(size=20,colour="black",family="Atlas Grotesk Web Bold",hjust=0.6),axis.title.x =element_blank())

r = r + theme(panel.grid.major=element_line(colour="white",size=0.5))
r = r + theme(panel.grid.minor=element_line(colour="white",size=0.5))
r = r + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=20),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))


give.n <- function(x){
  return(c(y = -0.2, label =length(x), size=8))
}


r = r+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))


ggsave(filename=paste(basename(args[1]),".domain.png",sep = ""),dpi=300)
