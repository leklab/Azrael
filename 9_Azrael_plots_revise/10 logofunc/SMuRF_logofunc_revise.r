#Rscript --vanilla SMuRF_logofunc_revise.r LARGE1_DiMSum_annotated_v7.tsv

library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)
library("ggrepel")

args = commandArgs(trailingOnly=TRUE)
df = read.csv(args[1], header=T, sep="\t")
df = df[df$confidence == "HIGH",]
df = df[df$LoGoFunc_GOF != "na",]


df$variantinfo <- paste("c.",df$WT_nt,df$site,df$Variant_nt, sep = "")


color_palette <- c("mediumspringgreen", "firebrick1", "lightblue")
df <- (df %>%
         arrange(factor(logofunc_pred, levels = c("Neutral", "LOF", "GOF"))))


h_line=0.1
# Main plot
scatterplot <- ggplot(df, aes(x = as.numeric(LoGoFunc_GOF), y = fitness_normalized, color=logofunc_pred,label=variantinfo)) +
  geom_point(size = 2) +geom_text_repel(aes(label=ifelse(as.numeric(df$LoGoFunc_GOF)>=0.3 & df$fitness_normalized>=h_line & df$logofunc_pred=="GOF",as.character(variantinfo),'')),size=3, face="bold")+
  geom_hline(yintercept = h_line, linetype = "dotted", color = "black",size=1) +
  geom_vline(xintercept = 0.3, linetype = "dotted", color = "black",size=1) + scale_color_manual(values = color_palette) +
  labs(x = "LoGoFunc_GOF score",
       y = "SMuRF score",
       color = "LoGoFunc prediction") +
  theme_minimal(base_family = "Helvetica") +
  theme(plot.title = element_text(size = 14, hjust = 0, vjust = 1, family = "Helvetica", face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = margin(5, 5, 5, 5)) + theme_bw()






variant_count_number=nrow(df)
scatterplot = scatterplot + labs(title="LARGE1", subtitle=paste(variant_count_number,"variants"))



scatterplot = scatterplot + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
scatterplot = scatterplot + theme(axis.text.x=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.x = element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))

scatterplot = scatterplot + theme(panel.grid.major=element_line(colour="white",size=0.5))
scatterplot = scatterplot + theme(panel.grid.minor=element_line(colour="white",size=0.5))

scatterplot = scatterplot + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))



scatterplot = scatterplot + theme(legend.text=element_text(colour="black",size=16))


scatterplot = scatterplot + theme(legend.position=c(0.25, 0.15),legend.title=element_text(colour="black",size=18))

scatterplot = scatterplot + guides(fill=guide_legend(title="Variant type"))

scatterplot = scatterplot + theme(legend.background=element_rect(linetype="solid",colour="black"))




# Print the scatterplot
ggsave(filename=paste(basename(args[1]),".style1_revise_hi_confi.png",sep = ""),dpi=300)
