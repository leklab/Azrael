#Rscript --vanilla SMuRF_phyloP_revise.r LARGE1_DiMSum_annotated_v5.tsv


library(ggplot2)



args = commandArgs(trailingOnly=TRUE)
df = read.csv(args[1],header=T,sep="\t")
df = df[!is.na(df$phyloP100way_vertebrate),]
df = df[df$classification == "Missense",]
df <- df[df$confidence == "HIGH",]

#1#



#1#
s = ggplot(df, aes(x=phyloP100way_vertebrate, y=fitness_normalized)) + theme_bw()





#1#
s = s + ylab("SMuRF score")+xlab("phyloP (100way vertebrate)")




s = s + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(colour="black"))
s = s + theme(axis.text.y=element_text(size=18,colour="black"),axis.title.y=element_text(size=18,colour="black",vjust=1.5))
s = s + theme(axis.text.x=element_text(size=18,colour="black"),axis.title.x= element_text(size=18,colour="black",vjust=1.5))


s = s + theme(panel.grid.major=element_line(colour="white",size=0.5))
s = s + theme(panel.grid.minor=element_line(colour="white",size=0.5))


variant_count_number=nrow(df)


########## FKRP or LARGE1 ######
#s = s + labs(title="FKRP Missense", subtitle=paste(variant_count_number,"variants"))
s = s + labs(title="LARGE1 Missense", subtitle=paste(variant_count_number,"variants"))



s = s + theme(axis.ticks=element_blank(),plot.title=element_text(size=14),plot.subtitle=element_text(size=11))
s = s  + geom_density_2d_filled(contour_var = "count")
s = s + geom_point(alpha=0.05, color="yellow")
s = s + theme(legend.text=element_text(colour="black",size=10))


s = s + theme(legend.title=element_text(colour="black",size=10))

s = s + guides(fill=guide_legend(title="Count"))

s = s + theme(legend.background=element_rect(linetype="solid",colour="black"))


#1
ggsave(filename=paste(basename(args[1]),"high_confidence.pdf",sep = ""))


#1#
cor.test(df$fitness_normalized, df$phyloP100way_vertebrate,method="spearman")
