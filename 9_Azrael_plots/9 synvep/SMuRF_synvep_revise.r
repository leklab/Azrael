#Rscript --vanilla SMuRF_synvep_revise.r FKRP_DiMSum_annotated_NEW_FINAL.tsv
#Rscript --vanilla SMuRF_synvep_revise.r LARGE1_DiMSum_annotated_v6.tsv

library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)
library("ggrepel")

args = commandArgs(trailingOnly=TRUE)
df = read.csv(args[1], header=T, sep="\t")
df = df[!is.na(df$synvep),]
df = df[df$confidence == "HIGH",]



synvep_cut=0.9
upc=1
botc=-2

df$variantinfo <- paste("c.",df$WT_nt,(as.numeric(df$site)),df$Variant_nt, sep = "")
df$color=ifelse(as.numeric(df$synvep)>=synvep_cut & df$fitness_normalized>=upc, "mediumspringgreen", "darkgreen")
df$color=ifelse(as.numeric(df$synvep)>=synvep_cut & df$fitness_normalized<=botc, "firebrick1", df$color)


# Main plot
scatterplot <- ggplot(df, aes(x = as.numeric(synvep), y = fitness_normalized,label=variantinfo, color=color)) +
  geom_point(size = 2) +geom_text_repel(aes(label=ifelse(as.numeric(df$synvep)>=synvep_cut & df$fitness_normalized>=upc,as.character(variantinfo),'')),size=3, fontface="bold", color=df$color)+ geom_text_repel(aes(label=ifelse(as.numeric(df$synvep)>=synvep_cut & df$fitness_normalized<=botc,as.character(variantinfo),'')),size=3, fontface="bold", color=df$color)+
geom_hline(yintercept = upc, linetype = "dotted", color = "black",size=1)+
  geom_vline(xintercept = synvep_cut, linetype = "dotted", color = "black",size=1)  +
geom_hline(yintercept = botc, linetype = "dotted", color = "black",size=1) + scale_color_manual(values = c("darkgreen", "firebrick1","mediumspringgreen")) +
  labs(x = "synVep score",
       y = "SMuRF score") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, hjust = 0, vjust = 1, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = margin(5, 5, 5, 5)) + theme_bw() + theme(legend.position = "none")






variant_count_number=nrow(df)
#gene#scatterplot = scatterplot + labs(title="FKRP", subtitle=paste(variant_count_number,"variants"))
scatterplot = scatterplot + labs(title="LARGE1", subtitle=paste(variant_count_number,"variants"))



scatterplot = scatterplot + theme(axis.text.y=element_text(size=18,colour="black"),axis.title.y=element_text(size=18,colour="black",vjust=1.5))
scatterplot = scatterplot + theme(axis.text.x=element_text(size=18,colour="black"),axis.title.x = element_text(size=18,colour="black",vjust=1.5))

scatterplot = scatterplot + theme(panel.grid.major=element_line(colour="white",size=0.5))
scatterplot = scatterplot + theme(panel.grid.minor=element_line(colour="white",size=0.5))

scatterplot = scatterplot + theme(axis.ticks=element_blank(),plot.title=element_text(size=14),plot.subtitle=element_text(size=11))





# Print the scatterplot
ggsave(filename=paste(basename(args[1]),".style1_revise_high_confi.pdf",sep = ""))
