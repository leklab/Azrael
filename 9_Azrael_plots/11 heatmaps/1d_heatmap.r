# 1d_heatmap.r
#
# This script processes the SMuRF variant count file to generate a heatmap of codon fitness scores.
# Usage: Rscript 1d_heatmap.r --input <input_file> --output <output_file>
#

library(ggplot2)
library(dplyr)
library(argparse)

parse_args <- function() {
  parser <- ArgumentParser(description = 'Generate heatmap from SMuRF variant count file.')
  parser$add_argument('-i', '--input', required = TRUE, help = 'Input SMuRF variant count file')
  parser$add_argument('-o', '--output', required = TRUE, help = 'Output PNG')
  return(parser$parse_args())
}

args <- parse_args()

df <- read.delim(args$input)

# Filter data
df <- df[df$classification == 'Missense',]
df <- df[df$confidence == "HIGH",]

# Aggregate fitness scores
new_df <- df %>%
  group_by(codon) %>%
  summarize(aggregated_fitness = mean(fitness_normalized, na.rm = TRUE))

new_df <- new_df %>%
  mutate(aggregated_fitness = pmax(pmin(aggregated_fitness, 2), -2))

# Generate heatmap
plot <- ggplot(new_df, aes(x = as.factor(codon), y = 1, fill = aggregated_fitness)) +
  geom_tile() +
  scale_fill_gradient2(low = "#B30000", mid = "#9dfa96", high = "#006400", midpoint = 0, limits = c(-2, 2),                      
                       name = "Mean missense score") +
  labs(x = "Codon", y = "") +
  scale_x_discrete(breaks = seq(0, max(new_df$codon), by = 100)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1/10,
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, 'cm'),
        panel.border = element_rect(fill = NA, colour = 'black', size = 0.5))

# Save plot
ggsave(args$output, plot = plot, device = "png", width = 8, height = 6, dpi = 600, bg = 'white')
