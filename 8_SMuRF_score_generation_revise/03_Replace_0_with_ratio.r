# 03_Replace_0_with_ratio.r
#
# This script processes the output from 02_Block_Parse.r and replaces values with no counts (0) with the ratio of the corresponding WT high to low.
# Usage: Rscript 03_Replace_0_with_ratio.r --prefix large1
# Output from this script will be in a format that is ready for DiMSum

library(argparse)

# Function to replace 0 values with the WT ratio
process_df <- function(df) {
  for (i in 2:nrow(df)) {
    for (j in seq(2, ncol(df), by = 2)) {
      high_col <- j     
      low_col <- j + 1  
      
      wt_high <- as.numeric(df[1, high_col])
      wt_low <- as.numeric(df[1, low_col])
      
      if (df[i, high_col] == 0) {
        df[i, high_col] <- round(wt_high / (wt_low + 1))
      }
      
      if (df[i, low_col] == 0) {
        df[i, low_col] <- 1
      }
    }
  }
  return(df)
}

# Argument parsing
parser <- ArgumentParser(description = 'Replace 0 values with the ratio of WT high to low.')
parser$add_argument('-p', '--prefix', type = 'character', required = TRUE, help = 'Prefix for the files from 02_Block_Parse.r')
args <- parser$parse_args()

# Generate file names based on the prefix
file_names <- list.files(pattern = paste0(args$prefix, "_variantCounts_block\\d+\\.txt"))

# Process each file
for (file_name in file_names) {
  df <- read.table(file_name, header = TRUE, sep = "\t")
  df <- process_df(df)
  write.table(df, file_name, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
