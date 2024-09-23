#Rscript --vanilla predict_heat_revise.r FKRP_DiMSum_annotated.tsv
#Rscript --vanilla predict_heat_revise.r LARGE1_DiMSum_annotated_v6.tsv

library(ggplot2)
library(reshape2)
library(dplyr)

# read file

args = commandArgs(trailingOnly=TRUE)
data = read.csv(args[1], header=T, sep="\t")
data = data[data$confidence == "HIGH",]

# only keep missense
data_missense = subset(data, classification == 'Missense')

# the 0s in these 2 columns should be NA

data_missense$EVE <- ifelse(data_missense$EVE == 0, NA, data_missense$EVE)
data_missense$mvp_rankscore <- ifelse(data_missense$mvp_rankscore == 0, NA, data_missense$mvp_rankscore)


# keep only the predictors and SMuRF

columns_to_keep <- c("fitness_normalized", "alphamissense", "primateAI_3D", "EVE", "revel_score", "MutScore","ESM1b_score","cadd","Maverick_RecScore", "meta_svm_score", "mvp_rankscore")
selected_data <- data_missense[, columns_to_keep]




# remove the NA rows

selected_data <- na.omit(selected_data)


# do as.numeric for all data

selected_data <- as.data.frame(sapply(selected_data, function(x) as.numeric(x)))



# Calculate the Spearman's rank correlation matrix
# Take the absolute value of the correlation matrix

cor_matrix <- cor(selected_data, method = "spearman")
abs_cor_matrix <- abs(cor_matrix)

colnames(abs_cor_matrix) <- c("SMuRF", "AlphaMissense", "PrimateAI-3D", "EVE", "REVEL", "MutScore", "ESM1b", "CADD","Maverick", "metaSVM", "MVP")
rownames(abs_cor_matrix) <- c("SMuRF", "AlphaMissense", "PrimateAI-3D", "EVE", "REVEL", "MutScore", "ESM1b", "CADD","Maverick", "metaSVM", "MVP")



# Reshape the correlation matrix for plotting
melted_data <- melt(abs_cor_matrix)



#remember to change gene name#
#remember to change gene name#
#remember to change gene name#
#remember to change gene name#
#remember to change gene name#

heatmap_plot <- ggplot(data = melted_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +  # You can customize colors here
  labs(x =element_blank(), y =element_blank(), fill = "Abs. Spearman's Rho", title="LARGE1") +theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(panel.border = element_rect(color = "black", size = 1))


heatmap_plot = heatmap_plot + theme(axis.text.y=element_text(size=15,colour="black",family="Atlas Grotesk Web Bold")) + theme(plot.title = element_text(size = 15, hjust = 0, vjust = 1, family = "Helvetica", face = "bold"),plot.margin = margin(5, 5, 5, 5))
heatmap_plot = heatmap_plot + theme(axis.text.x=element_text(size=15,colour="black",family="Atlas Grotesk Web Bold"))
heatmap_plot = heatmap_plot + theme(panel.grid.major=element_line(colour="white",size=0.5))
heatmap_plot = heatmap_plot + theme(panel.grid.minor=element_line(colour="white",size=0.5))

#heatmap_plot = heatmap_plot + theme(legend.text=element_text(colour="black",size=10))
#heatmap_plot = heatmap_plot + theme(legend.title=element_blank())
#heatmap_plot = heatmap_plot + theme(legend.background=element_blank())+ theme(legend.position="bottom")

heatmap_plot = heatmap_plot + guides(fill="none")


ggsave(filename=paste(basename(args[1]),".heatmap.png",sep = ""),dpi=300,height=12)
