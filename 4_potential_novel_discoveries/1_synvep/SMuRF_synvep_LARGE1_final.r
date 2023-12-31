#Rscript --vanilla SMuRF_synvep_LARGE1_final.r gnomad_large1_annotated_v3.tsv

library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)
library("ggrepel")

args = commandArgs(trailingOnly=TRUE)
df = read.csv(args[1], header=T, sep="\t")
df = df[!is.na(df$synvep),]



df$variantinfo <- paste("c.",df$WT_nt,df$site,df$Variant_nt, sep = "")
df$color=ifelse(as.numeric(df$synvep)>=0.9 & log2(df$functional_score)>=2, "mediumspringgreen", "darkgreen")
df$color=ifelse(as.numeric(df$synvep)>=0.9 & log2(df$functional_score)<=-2, "firebrick1", df$color)


# Main plot
scatterplot <- ggplot(df, aes(x = as.numeric(synvep), y = log2(functional_score),label=variantinfo, color=color)) +
  geom_point(size = 2) +geom_text_repel(aes(label=ifelse(as.numeric(df$synvep)>=0.9 & log2(df$functional_score)>=2,as.character(variantinfo),'')),size=3, fontface="bold", color=df$color, max.overlaps = Inf)+ geom_text_repel(aes(label=ifelse(as.numeric(df$synvep)>=0.9 & log2(df$functional_score)<=-2,as.character(variantinfo),'')),size=3, fontface="bold", color=df$color, max.overlaps = Inf)+
  geom_hline(yintercept = 2, linetype = "dotted", color = "black",size=1) +
  geom_vline(xintercept = 0.9, linetype = "dotted", color = "black",size=1)  +
geom_hline(yintercept = -2, linetype = "dotted", color = "black",size=1) + scale_color_manual(values = c("darkgreen", "firebrick1","mediumspringgreen")) +
  labs(x = "SynVep score",
       y = "log2(Functional score)") +
  theme_minimal(base_family = "Helvetica") +
  theme(plot.title = element_text(size = 14, hjust = 0, vjust = 1, family = "Helvetica", face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = margin(5, 5, 5, 5)) + theme_bw() + theme(legend.position = "none") + geom_text_repel(aes(label=ifelse(variantinfo=="c.T639C" | variantinfo=="c.A642G",as.character(variantinfo),'')),size=3, fontface="bold", color=df$color, max.overlaps = Inf)






variant_count_number=nrow(df)
scatterplot = scatterplot + labs(title="LARGE1", subtitle=paste(variant_count_number,"variants"))



scatterplot = scatterplot + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
scatterplot = scatterplot + theme(axis.text.x=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.x = element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))

scatterplot = scatterplot + theme(panel.grid.major=element_line(colour="white",size=0.5))
scatterplot = scatterplot + theme(panel.grid.minor=element_line(colour="white",size=0.5))

scatterplot = scatterplot + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))







# Print the scatterplot
ggsave(filename=paste(basename(args[1]),".style1.png",sep = ""),dpi=300)
