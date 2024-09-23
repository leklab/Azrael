# 05_format.r
#
# This script formats the output from DiMSum and performs confidence classification and normalization. The input prefix from DiMSum should be fitness_singles_b<block number>.txt.
# Usage: Rscript 05_format.r --blocksize <block_size> --gargamel <gargamel_file> --output <output_file>

library(dplyr)
library(argparse)

# Argument parsing
parser <- ArgumentParser(description = 'Format the output from DiMSum, perform confidence classification and normalization to get SMuRF scores')
parser$add_argument('-b', '--blocksize', type = 'integer', required = TRUE, help = 'Block size')
parser$add_argument('-g', '--gargamel', type = 'character', required = TRUE, help = 'Variant count file from Gargamel pipeline')
parser$add_argument('-o', '--output', type = 'character', required = TRUE, help = 'Output file name')

args <- parser$parse_args()

# List and order files based on the number of blocks
file_pattern <- "fitness_singles_b\\d+\\.txt"
files <- list.files(pattern = file_pattern)
files <- files[order(as.numeric(gsub("fitness_singles_b(\\d+)\\.txt", "\\1", files)))]

dfs <- list()

for (file in files) {
  dfs[[file]] <- read.table(file, header = TRUE)
}

for (i in 2:length(dfs)) {
  increment <- (i - 1) * args$blocksize
  dfs[[i]]$Pos <- dfs[[i]]$Pos + increment
}

df <- do.call(rbind, dfs)
rownames(df) <- NULL
df$WT_nt <- toupper(df$WT_nt)
df$Mut <- toupper(df$Mut)
df$variant_key <- paste(df$Pos, df$WT_nt, df$Mut, sep = ":")

# Read original Gargamel variant file
original <- read.delim(args$gargamel, header = TRUE)
original$variant_key <- paste(original$site, original$WT_nt, original$Variant_nt, sep = ":")

df <- df[order(df$variant_key), ]
original <- original[order(original$variant_key), ]
merged_df <- merge(df, original, by = "variant_key", all.x = TRUE)

df$classification <- merged_df$classification
df <- df[order(df$Pos), ]

# Confidence classification
Q1 <- quantile(df$sigma, 0.25)
Q3 <- quantile(df$sigma, 0.75)
IQR <- Q3 - Q1
upper_bound <- Q3 + 1.5 * IQR
df$confidence <- ifelse(df$sigma <= upper_bound, "HIGH", "LOW")

# Merge original formatting from Gargamel pipeline
df$nt_seq <- NULL
df$Nham_nt <- NULL

merged_df <- merge(df, original[, c("chr", "hg38_site", "site", "codon", "Variant_nt", "WT_AA", "Variant_AA", "variant_key", 'Clinvar_clinical_significance')], by = "variant_key", all.x = TRUE)

merged_df$nt_seq <- NULL
merged_df$Pos <- NULL
merged_df$Mut <- NULL
merged_df$variant_key <- NULL

# Reorder columns
cols_to_reorder <- c("chr", "hg38_site", "site", "codon", "WT_nt", "Variant_nt", "WT_AA", "Variant_AA", "Clinvar_clinical_significance")
merged_df <- merged_df[, c(cols_to_reorder, setdiff(names(merged_df), cols_to_reorder))]

merged_df <- merged_df[order(merged_df$Variant_nt), ]
merged_df <- merged_df[order(merged_df$site), ]

# Normalize scores based on synonymous
synonymous_subset <- merged_df %>% filter(classification == "Synonymous")

synonymous_subset <- synonymous_subset %>%
  mutate(fitness_normal_scale = exp(fitness))

mean_fitness_by_block <- synonymous_subset %>%
  group_by(Block) %>%
  summarise(mean_fitness_normal = mean(fitness_normal_scale, na.rm = TRUE))

merged_df <- merged_df %>%
  left_join(mean_fitness_by_block, by = "Block") %>%
  rename(scaling_factor = mean_fitness_normal)

merged_df <- merged_df %>%
  mutate(fitness_normalized = log(exp(fitness) / scaling_factor))

merged_df$fitness <- log2(exp(merged_df$fitness))
merged_df$fitness_normalized <- log2(exp(merged_df$fitness_normalized))

# Write the final output
write.table(merged_df, file = args$output, row.names = FALSE, sep = "\t", quote = FALSE)

