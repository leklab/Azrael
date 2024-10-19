# 2d_heatmap.r
#
# This script processes the SMuRF variant count file to generate a 2D heatmap.
# Usage: Rscript 2d_heatmap.r --input <input_file> --codon_segment <codon_segment_size> --output <output_prefix>
#

library(ggplot2)
library(dplyr)
library(cowplot)
library(svglite)
library(argparse)

parse_args <- function() {
  parser <- ArgumentParser(description = 'Generate 2D heatmap from SMuRF variant count file.')
  parser$add_argument('-i', '--input', required = TRUE, help = 'Input SMuRF variant count file')
  parser$add_argument('-c', '--codon_segment', type = "integer", default = 100, help = 'Size of codon segment per row in heatmap (default: 100)')
  parser$add_argument('-o', '--output', required = TRUE, help = 'Output file prefix for the PNG and SVG')
  return(parser$parse_args())
}

args <- parse_args()

df <- read.delim(args$input)

# Filter data
df_high_conf <- df %>%
  filter(classification == 'Missense', confidence == "HIGH")

df_low_conf <- df %>%
  filter(classification == 'Missense', confidence == "LOW")

# Define amino acid groups
nonpolar_aliphatic <- c('G', 'A', 'V', 'L', 'M', 'I')
polar_uncharged <- c('S', 'T', 'C', 'P', 'N', 'Q')
positively_charged <- c('K', 'R', 'H')
negatively_charged <- c('D', 'E')
nonpolar_aromatic <- c('F', 'Y', 'W')

df_high_conf <- df_high_conf %>%
  mutate(AA_Group = case_when(
    Variant_AA %in% nonpolar_aliphatic ~ 'Nonpolar, aliphatic',
    Variant_AA %in% polar_uncharged ~ 'Polar, uncharged',
    Variant_AA %in% positively_charged ~ 'Positively charged',
    Variant_AA %in% negatively_charged ~ 'Negatively charged',
    Variant_AA %in% nonpolar_aromatic ~ 'Nonpolar, aromatic',
    TRUE ~ NA_character_
  ))

df_high_conf <- df_high_conf %>%
  mutate(codon = as.factor(codon))

summary_df <- df_high_conf %>%
  group_by(codon, AA_Group) %>%
  summarise(mean_fitness = mean(fitness_normalized, na.rm = TRUE))

summary_df <- summary_df %>%
  mutate(mean_fitness = pmax(pmin(mean_fitness, 2), -2))

wt_df <- df_high_conf %>%
  distinct(codon, WT_AA) %>%
  mutate(AA_Group_wt = case_when(
    WT_AA %in% nonpolar_aliphatic ~ 'Nonpolar, aliphatic',
    WT_AA %in% polar_uncharged ~ 'Polar, uncharged',
    WT_AA %in% positively_charged ~ 'Positively charged',
    WT_AA %in% negatively_charged ~ 'Negatively charged',
    WT_AA %in% nonpolar_aromatic ~ 'Nonpolar, aromatic',
    TRUE ~ NA_character_
  ))

template_df <- expand.grid(codon = unique(summary_df$codon), AA_Group = unique(summary_df$AA_Group))

summary_df <- merge(summary_df, template_df, by = c("codon", "AA_Group"), all = TRUE)

summary_df <- merge(summary_df, wt_df, by = "codon", all.x = TRUE)
summary_df$wt_match <- ifelse(summary_df$AA_Group == summary_df$AA_Group_wt, TRUE, FALSE)

summary_df$low_confidence <- FALSE

low_conf_mark <- df_low_conf %>%
  mutate(AA_Group = case_when(
    Variant_AA %in% nonpolar_aliphatic ~ 'Nonpolar, aliphatic',
    Variant_AA %in% polar_uncharged ~ 'Polar, uncharged',
    Variant_AA %in% positively_charged ~ 'Positively charged',
    Variant_AA %in% negatively_charged ~ 'Negatively charged',
    Variant_AA %in% nonpolar_aromatic ~ 'Nonpolar, aromatic',
    TRUE ~ NA_character_
  )) %>%
  distinct(codon, AA_Group) %>%
  mutate(low_confidence = TRUE)

summary_df$codon <- as.character(summary_df$codon)
low_conf_mark$codon <- as.character(low_conf_mark$codon)

summary_df <- merge(summary_df, low_conf_mark, by = c("codon", "AA_Group"), all.x = TRUE, suffixes = c("", "_low_conf"))

summary_df$low_confidence <- ifelse(is.na(summary_df$low_confidence_low_conf), FALSE, summary_df$low_confidence_low_conf)

summary_df <- summary_df %>%
  select(-low_confidence_low_conf)

summary_df$codon <- factor(summary_df$codon, levels = sort(unique(as.numeric(as.character(summary_df$codon)))))

order_groups <- c("Negatively charged", "Positively charged", "Polar, uncharged", "Nonpolar, aromatic", "Nonpolar, aliphatic")

summary_df$AA_Group <- factor(summary_df$AA_Group, levels = order_groups)

plot_segment <- function(data, start, end) {
  ggplot(data, aes(x = codon, y = AA_Group, fill = mean_fitness)) +
    geom_tile(color = "white") +
    geom_point(aes(size = ifelse(wt_match, "wt_dot", "no_dot")), shape = 16, color = "black") + 
    geom_point(data = subset(data, is.na(mean_fitness) & low_confidence), 
               aes(x = codon, y = AA_Group), shape = 4, size = 3, color = "red") +
    scale_size_manual(values = c(wt_dot = 1, no_dot = NA), guide = "none") +
    coord_equal() +
    scale_fill_gradient2(low = "#B30000", mid = "#9dfa96", high = "#006400", midpoint = 0, limits = c(-2, 2), 
                         name = "Mean missense score") +
    scale_x_discrete(breaks = seq(start, end, by = 20)) +
    scale_y_discrete(expand = c(0, 0), limits = rev(order_groups)) +
    theme_classic() +
    labs(x = "",
         y = "") +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          axis.title.x = element_text(size = 16),
          plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm"),
          legend.position = "none")
}

# Define the number of codons per segment
codons_per_segment <- args$codon_segment

# Determine the number of codons and segments
num_codons <- max(as.numeric(as.character(summary_df$codon)))
num_segments <- ceiling(num_codons / codons_per_segment)

# Generate and combine plots for each segment
plots <- list()
for (i in 1:num_segments) {
  start <- (i - 1) * codons_per_segment + 1
  end <- min(i * codons_per_segment, num_codons)
  segment_df <- summary_df %>%
    filter(as.numeric(as.character(codon)) >= start & as.numeric(as.character(codon)) <= end)
  plots[[i]] <- plot_segment(segment_df, start, end)
}

# Combine all the plots and add a single legend at the bottom
combined_plot <- plot_grid(plotlist = plots, ncol = 1)

# Extract legend from one of the plots
legend <- get_legend(
  ggplot(summary_df, aes(x = codon, y = AA_Group, fill = mean_fitness)) +
    geom_tile(color = "white") +
    coord_equal() +
    scale_fill_gradient2(low = "#B30000", mid = "#9dfa96", high = "#006400", midpoint = 0, limits = c(-2, 2), 
                         name = "Mean missense score") +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          legend.key.size = unit(1, "lines"), 
          legend.spacing.x = unit(0.5, "lines"))
)

# Combine the plots and the legend
final_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(12, 1))

# Save the final plot with a white background
output_png <- paste0(args$output, ".png")
output_svg <- paste0(args$output, ".svg")

ggsave(output_png, plot = final_plot, device = "png", width = 20, height = 10, dpi = 600, bg = 'white')
ggsave(output_svg, plot = final_plot, device = "svg", width = 20, height = 10, bg = "white")
