# Rscript --vanilla SMuRF_clinvar_large1_revise_main.r LARGE1_DiMSum_annotated_v2.tsv



library(ggplot2)
library(ggsignif)

#wilcox.test.default()
args = commandArgs(trailingOnly=TRUE)


all_blocks = read.csv(args[1],header=T,sep="\t")
all_blocks <- all_blocks[all_blocks$confidence == "HIGH",]

all_blocks_missense = subset(all_blocks, classification == 'Missense')


all_blocks_clinvar = all_blocks_missense[!is.na(all_blocks_missense$Clinvar_clinical_significance) & all_blocks_missense$Clinvar_clinical_significance!='Conflicting interpretations of pathogenicity' & all_blocks_missense$Clinvar_clinical_significance!='no interpretation for the single variant',]

all_blocks_clinvar$simplified_labels="VUS"
all_blocks_clinvar$simplified_labels=ifelse(all_blocks_clinvar$Clinvar_clinical_significance=='Benign'| all_blocks_clinvar$Clinvar_clinical_significance=='Benign/Likely benign'| all_blocks_clinvar$Clinvar_clinical_significance=='Likely benign', "B/LB", all_blocks_clinvar$simplified_labels)
all_blocks_clinvar$simplified_labels=ifelse(all_blocks_clinvar$Clinvar_clinical_significance=='Pathogenic'| all_blocks_clinvar$Clinvar_clinical_significance=='Pathogenic/Likely pathogenic'| all_blocks_clinvar$Clinvar_clinical_significance=='Likely pathogenic', "P/LP", all_blocks_clinvar$simplified_labels)



level_order <- c('B/LB','VUS','P/LP')


p = ggplot(all_blocks_clinvar, aes(x=factor(simplified_labels,level_order), y=fitness_normalized,color=simplified_labels)) + theme_bw()
p = p+scale_color_manual(values=c("darkgreen","red","lightblue"))
p = p+ theme(legend.position="none")



p = p + geom_violin(lwd=0.6)+geom_signif(comparisons = list(c("B/LB", "VUS")),y_position=2.5,size=1,textsize=5,color="black")+geom_signif(comparisons = list(c("B/LB", "P/LP")),y_position=3,size=1,textsize=5,color="black")+geom_signif(comparisons = list(c("VUS", "P/LP")),y_position=2,size=1,textsize=5,color="black")
p = p + geom_boxplot(width=0.1,outlier.size=3)
p = p + scale_y_continuous("SMuRF score") + theme(axis.text.x=element_text(size=12), plot.margin = unit(c(0.5, 1, 1, 1), "cm"))
p = p + geom_hline(yintercept=0, linetype="dashed", color="black",size=0.5) + scale_x_discrete(name="")


give.n <- function(x){
  return(c(y = -4.8, label =length(x), size=10))
}



p=p+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))

p = p + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(color="black",size=1.5))
p = p + theme(axis.text.y=element_text(size=20,colour="black",family="Atlas Grotesk Web Bold"),axis.title.y=element_text(size=25,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))


p = p + scale_x_discrete(labels= level_order)
p = p + theme(axis.text.x=element_text(size=25,colour="black",family="Atlas Grotesk Web Bold"),axis.title.x = element_blank())



p = p + theme(panel.grid.major=element_line(colour="white",size=0.5))
p = p + theme(panel.grid.minor=element_line(colour="white",size=0.5))
p = p + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=20),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=18))

variant_count_number=nrow(all_blocks_clinvar)
p = p + labs(title="LARGE1 Missense", subtitle=paste(variant_count_number,"variants"))



ggsave(filename="SMuRF_clinvar_large1_revise_main_high_confidence.png",dpi=300)
