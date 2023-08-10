#Rscript --vanilla SMuRF_gnomad_LARGE1_final.r gnomad_large1_annotated_v2.tsv
#DOUBLE CHECK TO MAKE SURE LEGEND DOES NOT MASK ANY DATA POINT

library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
df = read.csv(args[1], header=T, sep="\t")


df <- df[df$gnomad_v3_af != 0,]

min_allele_freq <- min(df$gnomad_v3_af)
max_allele_freq <- max(df$gnomad_v3_af)

color_palette <- c("lightblue", "red","darkred","darkgreen")

# arrange

df <- (df %>%
         arrange(factor(classification, levels = c("Missense", "Synonymous", "Nonsense","Start-loss"))))







# Main plot
scatterplot <- ggplot(df, aes(x = gnomad_v3_af, y = log2(functional_score), color = classification)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black",size=1) +
  scale_x_continuous(trans = "log10", limits = c(min_allele_freq, max_allele_freq),
                     labels = function(x) parse(text = format(x, scientific = FALSE))) +  #Log transform the x-axis scale (but not the x-axis values)
  scale_color_manual(values = color_palette) + # Apply the color palette
  labs(x = "gnomAD v3 Allele Frequency",
       y = "log2(Functional score)",
       color = "Classification") +
  theme_minimal(base_family = "Helvetica") +
  theme(plot.title = element_text(size = 14, hjust = 0, vjust = 1, family = "Helvetica", face = "bold"), # Adjustments of plot title font size, alignment, and font style
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = margin(5, 5, 5, 5)) + theme_bw()


#Add style

variant_count_number=nrow(df)
scatterplot = scatterplot + labs(title="LARGE1", subtitle=paste(variant_count_number,"variants"))



scatterplot = scatterplot + theme(axis.text.y=element_text(size=18,colour="black",family="Decima Mono Pro"),axis.title.y=element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))
scatterplot = scatterplot + theme(axis.text.x=element_text(size=18,colour="black"),axis.title.x = element_text(size=18,colour="black",family="Atlas Grotesk Web Bold",vjust=1.5))

scatterplot = scatterplot + theme(panel.grid.major=element_line(colour="white",size=0.5))
scatterplot = scatterplot + theme(panel.grid.minor=element_line(colour="white",size=0.5))

scatterplot = scatterplot + theme(axis.ticks=element_blank(),plot.title=element_text(family="Atlas Grotesk Web Bold", size=14),plot.subtitle=element_text(family="Atlas Grotesk Web Light", size=11))
scatterplot = scatterplot + theme(legend.text=element_text(colour="black",size=16))


scatterplot = scatterplot + theme(legend.position=c(0.825, 0.15),legend.title=element_text(colour="black",size=18))

scatterplot = scatterplot + guides(fill=guide_legend(title="Variant type"))

scatterplot = scatterplot + theme(legend.background=element_rect(linetype="solid",colour="black"))










# Add a grey rectangle box to the main plot to indicate low allele counts
scatterplot <- scatterplot +
  annotate("rect", xmin = 0, xmax = 1.5e-05, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5)


# Print the scatterplot
ggsave(filename=paste(basename(args[1]),".gnomad.png",sep = ""),dpi=300)
