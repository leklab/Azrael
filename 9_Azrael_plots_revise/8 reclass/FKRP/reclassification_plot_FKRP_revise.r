#Rscript --vanilla reclassification_plot_FKRP_revise.r count_FKRP_revise

library(ggplot2)
library(ggalluvial)


args = commandArgs(trailingOnly=TRUE)
data = read.csv(args[1], header=T, sep="\t")
data$Counts=as.numeric(data$Counts)

#1#
data=data[data$Clinvar_clinical_significance!="Unclassified", ]

#2#
#2#data=data[data$Clinvar_clinical_significance=="Unclassified", ]

level_order <- c('B/LB','VUS','P/LP', 'Unclassified')
level_order2 <- c('SMuRF Benign','SMuRF Mild', "SMuRF Intermediate", "SMuRF Severe")





p = ggplot(data, aes(axis2 = factor(SMuRF_reclassification,level_order2),axis1=factor(Clinvar_clinical_significance,level_order), y = Counts)) + theme_void()


p = p + geom_alluvium(aes(fill=SMuRF_reclassification),width=1/10) + geom_stratum(aes(fill=SMuRF_reclassification),width=1/5)

p = p + theme(plot.title=element_text(family="Atlas Grotesk Web Bold", size=28,hjust=0.5))


p = p + geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=6,fontface="bold", color="white")
p = p + scale_x_discrete(limits = c("Clinvar_clinical_significance", "SMuRF_reclassification"), expand = c(0, 0))
p = p + theme(legend.position = "none") + scale_fill_manual(values=c("darkblue","violetred3","pink","darkred"))

p = p + theme(axis.text.y=element_text(size=20,colour="black",family="Decima Mono Pro",face="bold"))


#1#
ggsave(filename=paste(basename(args[1]),".reclass_revise.png",sep = ""),bg='white',dpi=300,width=16,height=8)

#2#
#2#ggsave(filename=paste(basename(args[1]),".class_revise.png",sep = ""),bg='white',dpi=300,width=16,height=8)
