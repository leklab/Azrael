# 03_Replace_0_with_WT.r
#
# This script processes the output from 02_Block_Parse.r (remember to add in the wild-type values!) and replaces values with no counts (0) with the WT value.
# Usage: Rscript 03_Replace_0_with_WT.r --prefix large1
# Output from this script will be in a format that is ready for DiMSum


library(argparse)

# Function to replace 0 values with the WT values (first row)
process_df <- function(df) {
  for (col in names(df)[-1]) {
    df[[col]][df[[col]] == 0] <- df[[col]][1]
  }
  return(df)
}

count_rows_with_zero <- function(df) {
  row_has_zero <- apply(df[-1], 1, function(row) any(row == 0))
  return(sum(row_has_zero))
}

# Argument parsing
parser <- ArgumentParser(description = 'Replace values with no counts (0) with the WT value.')
parser$add_argument('-p', '--prefix', type = 'character', required = TRUE, help = 'Prefix for the files from 02_Block_Parse.r')
args <- parser$parse_args()

# Generate file names based on the prefix
file_names <- list.files(pattern = paste0(args$prefix, "_variantCounts_block\\d+\\.txt"))

# Process each file
file_count <- length(file_names)
for (file_name in file_names) {
  df <- read.table(file_name, header = TRUE, sep = "\t")
  df <- process_df(df)
  write.table(df, file_name, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
