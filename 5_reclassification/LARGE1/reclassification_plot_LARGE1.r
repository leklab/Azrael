#Rscript --vanilla reclassification_plot_LARGE1.r count_LARGE1

library(ggplot2)
library(ggalluvial)




args = commandArgs(trailingOnly=TRUE)

data = read.csv(args[1], header=T, sep="\t")

data$Counts=as.numeric(data$Counts)

#1#
data=data[data$Clinvar_clinical_significance!="Unclassified", ]

#2#
#2#data=data[data$Clinvar_clinical_significance=="Unclassified", ]

level_order <- c('P/LP','B/LB','VUS', 'Unclassified')

level_order2 <- c('SMuRF Benign','SMuRF Pathogenic')

level_order3 <- c('High', 'Medium', 'Low')




p = ggplot(data, aes(axis3 = factor(confidence, level_order3), axis2 = factor(SMuRF_reclassification,level_order2),axis1=factor(Clinvar_clinical_significance,level_order), y = Counts)) + theme_void()


p = p + geom_alluvium(aes(fill=SMuRF_reclassification),width=1/8) + geom_stratum(aes(fill=SMuRF_reclassification),width=1/7)

p = p + theme(plot.title=element_text(family="Atlas Grotesk Web Bold", size=28,hjust=0.5))


p = p + geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=6,fontface="bold", color="white")
p = p + scale_x_discrete(limits = c("Clinvar_clinical_significance", "SMuRF_reclassification","confidence"), expand = c(0.08, 0.0))
p = p + theme(legend.position = "none") + scale_fill_manual(values=c("darkblue","darkred"))

p = p + theme(axis.text.y=element_text(size=20,colour="black",family="Decima Mono Pro",face="bold"))


#1#
ggsave(filename=paste(basename(args[1]),".reclass.png",sep = ""),bg='white',dpi=300,width=35,height=12)

#2#
#2#ggsave(filename=paste(basename(args[1]),".class.png",sep = ""),bg='white',dpi=300,width=35,height=12)
