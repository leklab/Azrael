# 02_Block_Parse.r
#
# This script processes the output from 01_Variant_Prep_<gene>.py and splits the nucleotide sequences into blocks for DiMSum
# Usage: Rscript 02_Block_Parse.r --variant large1_variantCounts_site.txt --wildtype large1_sequence.txt --block 227 --prefix large1
# Remember to add in the WT values in the first row for the output files for this script before using it for 03_Replace_0_with_WT.r

library(dplyr)
library(argparse)

# Argument parsing
parser <- ArgumentParser(description = 'Split nucleotide sequences from into blocks for DiMSum')
parser$add_argument('-v', '--variant', type = 'character', required = TRUE, help = 'Output from 01_Variant_Prep_<gene>.py')
parser$add_argument('-w', '--wildtype', type = 'character', required = TRUE, help = 'Wild-type nucleotide sequence')
parser$add_argument('-b', '--block', type = 'integer', required = TRUE, help = 'Block size')
parser$add_argument('-p', '--prefix', type = 'character', required = TRUE, help = 'Output file prefix')

args <- parser$parse_args()

# Read the variant file
df <- read.delim(args$variant)

# Read the WT sequence
wt <- readLines(args$wildtype)
wt <- paste(wt, collapse = "")

# Define the intervals and blocks
block_size <- args$block
intervals <- seq(1, nchar(wt), by = block_size)
df$block <- cut(df$site, breaks = seq(0, nchar(wt) + 1, by = block_size), labels = FALSE)

# Split the nt_seq into blocks
df <- df %>%
  rowwise() %>%
  mutate(nt_seq_block = substring(nt_seq, intervals[block], 
                                  ifelse(block < length(intervals), intervals[block + 1] - 1, nchar(nt_seq))),
         nt_seq_block_char = nchar(nt_seq_block)) %>%
  ungroup()

df_list <- split(df, df$block)

# Generate the wild-type blocks
blocks <- sapply(seq_along(intervals), function(i) {
  start <- intervals[i]
  end <- ifelse(i < length(intervals), intervals[i + 1] - 1, nchar(wt))
  substring(wt, start, end)
})

# Write the output files
for (i in seq_along(df_list)) {
  df_list[[i]] <- df_list[[i]] %>%
    select(-c(site, block, nt_seq_block_char, nt_seq)) %>%
    rename(nt_seq = nt_seq_block) %>%
    select(nt_seq, everything()) %>%
    add_row(nt_seq = blocks[i], .before = 1)
  
  output_file <- paste0(args$prefix, "_variantCounts_block", i, ".txt")
  write.table(df_list[[i]], file = output_file, sep = "\t", row.names = FALSE)
}
