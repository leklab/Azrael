#Rscript --vanilla SMuRF_gnomadv4_FKRP_revise_exogeno.r FKRP_DiMSum_annotated.tsv


library(ggplot2)
library(scales)
library(RColorBrewer)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
df = read.csv(args[1], header=T, sep="\t")


df <- df[df$gnomad_v4_genome_exome_af != 0,]
df <- df[df$confidence == "HIGH",]

min_allele_freq <- min(df$gnomad_v4_genome_exome_af)-0.06
max_allele_freq <- max(df$gnomad_v4_genome_exome_af)+0.06

color_palette <- c("lightblue", "red","darkred","darkgreen")

# arrange

df <- (df %>%
         arrange(factor(classification, levels = c("Missense", "Synonymous", "Nonsense","Start-loss"))))






# Main plot
scatterplot <- ggplot(df, aes(x = gnomad_v4_genome_exome_af, y = fitness_normalized, color = classification)) +
  geom_jitter(size = 3,width = 0.05, height = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",size=1) +
  scale_x_continuous(trans = "log10", limits = c(min_allele_freq, max_allele_freq),
                     labels = function(x) parse(text = format(x, scientific = FALSE))) +  #Log transform the x-axis scale (but not the x-axis values)
  scale_color_manual(values = color_palette) + # Apply the color palette
  labs(x = "gnomAD v4 Allele Frequency",
       y = "SMuRF score",
       color = "Classification") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, hjust = 0, vjust = 1, face = "bold"), # Adjustments of plot title font size, alignment, and font style
        plot.margin = margin(5, 5, 5, 5)) + theme_bw()


#Add style

variant_count_number=nrow(df)
scatterplot = scatterplot + labs(title="FKRP", subtitle=paste(variant_count_number,"variants"))
scatterplot = scatterplot + theme(panel.border = element_rect(color = "black", size=1.5))


scatterplot = scatterplot + theme(axis.text.y=element_text(size=20,colour="black"),axis.title.y=element_text(size=22,colour="black",vjust=1.5))
scatterplot = scatterplot + theme(axis.text.x=element_text(size=20,colour="black"),axis.title.x = element_text(size=22,colour="black",vjust=1.5))

scatterplot = scatterplot + theme(panel.grid.major=element_line(colour="white",size=0.5))
scatterplot = scatterplot + theme(panel.grid.minor=element_line(colour="white",size=0.5))

scatterplot = scatterplot + theme(axis.ticks=element_blank(),plot.title=element_text(size=20),plot.subtitle=element_text(size=18))
scatterplot = scatterplot + theme(legend.text=element_text(colour="black",size=16))


scatterplot = scatterplot + theme(legend.position=c(0.84, 0.15),legend.title=element_text(colour="black",size=18))

scatterplot = scatterplot + guides(fill=guide_legend(title="Variant type"))

scatterplot = scatterplot + theme(legend.background=element_rect(linetype="solid",colour="black", fill="NA"))










# Add a grey rectangle box to the main plot to indicate low allele counts
scatterplot <- scatterplot +
  annotate("rect", xmin = 0, xmax = 1.5e-05, ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.5)



# Print the scatterplot
ggsave(filename=paste(basename(args[1]),".gnomadv4_exogeno_jitter_high_confidence.pdf",sep = ""))
