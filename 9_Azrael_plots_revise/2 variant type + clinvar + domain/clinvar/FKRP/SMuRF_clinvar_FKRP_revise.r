# Rscript --vanilla SMuRF_clinvar_FKRP_revise.r FKRP_DiMSum_annotated.tsv




library(ggplot2)
args = commandArgs(trailingOnly=TRUE)


all_blocks = read.csv(args[1],header=T,sep="\t")
all_blocks <- all_blocks[all_blocks$confidence == "HIGH",]


all_blocks_clinvar = all_blocks[!is.na(all_blocks$Clinvar_clinical_significance) & all_blocks$Clinvar_clinical_significance!='Conflicting interpretations of pathogenicity' & all_blocks$Clinvar_clinical_significance!='no interpretation for the single variant',]


level_order <- c('Benign', 'Benign/Likely benign', 'Likely benign','Uncertain significance','Likely pathogenic','Pathogenic/Likely pathogenic','Pathogenic')


p = ggplot(all_blocks_clinvar, aes(x=factor(Clinvar_clinical_significance,level_order), y=fitness_normalized,color=Clinvar_clinical_significance)) + theme_bw()
p = p + scale_color_manual(values=c("darkgreen", "darkgreen","darkgreen","red","red","red","lightblue"))
p = p + theme(legend.position="none")





p = p + geom_violin(lwd=0.6)
p = p + geom_boxplot(width=0.1,outlier.size=3)
p = p + scale_y_continuous("SMuRF score") + theme(axis.text.x=element_text(size=12), plot.margin = unit(c(0.5, 1, 1, 1), "cm"))

p = p + geom_hline(yintercept=0, linetype="dashed", color="black",size=0.5) + scale_x_discrete(name="") + ggtitle("Clivar interpretation - All blocks")


give.n <- function(x){
  return(c(y = 3.4, label =length(x), size=10))
}



p=p+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))

p = p + theme(panel.background=element_rect(fill="white"), plot.background=element_rect(fill="white"),panel.border=element_rect(color="black",size=1.5))
p = p + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))

simplified_labels=c('B', 'B/LB', 'LB','VUS','LP', 'P/LP','P')
p = p + scale_x_discrete(labels= simplified_labels)
p = p + theme(axis.text.x=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold"),axis.title.x =element_blank())



p = p + theme(panel.grid.major=element_line(colour="white",size=0.5))
p = p + theme(panel.grid.minor=element_line(colour="white",size=0.5))
p = p + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=20),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=18))

p=p+stat_summary(fun.data = give.n, geom = "text", fun.y = median, position = position_dodge(width = 0.75))

variant_count_number=nrow(all_blocks_clinvar)
p = p + labs(title="FKRP", subtitle=paste(variant_count_number,"variants"))


ggsave(filename="SMuRF_clinvar_FKRP_revise_high_confidence.png",dpi=300, width=7)
