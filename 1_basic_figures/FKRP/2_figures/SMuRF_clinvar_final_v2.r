#Rscript --vanilla SMuRF_clinvar_final_v2.r normalized_FKRP_block1_func normalized_FKRP_block2_func  normalized_FKRP_block3_func normalized_FKRP_block4_func normalized_FKRP_block5_func normalized_FKRP_block6_func

library(ggplot2)
library(ggsignif)

#wilcox.test.default()

args = commandArgs(trailingOnly=TRUE)


block1 = read.csv(args[1],header=T,sep="\t")
block2 = read.csv(args[2],header=T,sep="\t")
block3 = read.csv(args[3],header=T,sep="\t")
block4 = read.csv(args[4],header=T,sep="\t")
block5 = read.csv(args[5],header=T,sep="\t")
block6 = read.csv(args[6],header=T,sep="\t")





block1$block = 'block1'
block2$block = 'block2'
block3$block = 'block3'
block4$block = 'block4'
block5$block = 'block5'
block6$block = 'block6'





all_blocks = rbind(block1,block2,block3,block4,block5,block6)
all_blocks_clinvar = all_blocks[!is.na(all_blocks$Clinvar_clinical_significance) & all_blocks$Clinvar_clinical_significance!='Conflicting interpretations of pathogenicity' & all_blocks$Clinvar_clinical_significance!='no interpretation for the single variant',]

all_blocks_clinvar$simplified_labels="VUS"
all_blocks_clinvar$simplified_labels=ifelse(all_blocks_clinvar$Clinvar_clinical_significance=='Benign'| all_blocks_clinvar$Clinvar_clinical_significance=='Benign/Likely benign'| all_blocks_clinvar$Clinvar_clinical_significance=='Likely benign', "B/LB", all_blocks_clinvar$simplified_labels)
all_blocks_clinvar$simplified_labels=ifelse(all_blocks_clinvar$Clinvar_clinical_significance=='Pathogenic'| all_blocks_clinvar$Clinvar_clinical_significance=='Pathogenic/Likely pathogenic'| all_blocks_clinvar$Clinvar_clinical_significance=='Likely pathogenic', "P/LP", all_blocks_clinvar$simplified_labels)


level_order <- c('B/LB','VUS','P/LP')


p = ggplot(all_blocks_clinvar, aes(x=factor(simplified_labels,level_order), y=log2(functional_score),fill=simplified_labels)) + theme_bw()
p = p+scale_fill_manual(values=c("darkgreen","red","lightblue"))
p = p+ theme(legend.position="none")+geom_signif(comparisons = list(c("B/LB", "VUS")),y_position=4)+geom_signif(comparisons = list(c("B/LB", "P/LP")),y_position=4.4)+geom_signif(comparisons = list(c("VUS", "P/LP")),y_position=3.6)




p = p + geom_boxplot()
p = p + scale_y_continuous("log2(Functional score)") + theme(axis.text.x=element_text( size=12), plot.margin = unit(c(0.5, 1, 1, 1), "cm"))
p = p + geom_hline(yintercept=0, linetype="dashed", color="blue") + scale_x_discrete(name="") + ggtitle("Clivar interpretation - All blocks")


give.n <- function(x){
  return(c(y = -4.8, label =length(x), size=5))
}



p=p+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))

p = p + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_blank())
p = p + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))


p = p + scale_x_discrete(labels= level_order)
p = p + theme(axis.text.x=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.x = element_blank())



p = p + theme(panel.grid.major=element_line(colour="white",size=0.5))
p = p + theme(panel.grid.minor=element_line(colour="white",size=0.5))
p = p + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))

variant_count_number=nrow(all_blocks_clinvar)
p = p + labs(title="FKRP", subtitle=paste(variant_count_number,"variants"))



ggsave(filename="SMuRF_clinvar_v2.png",dpi=300)
