#!/usr/bin/env Rscript

#***** DOES NOT WORK *****#


# Load required packages
suppressPackageStartupMessages({
  library(pdftools)
  library(stringr)
  library(dplyr)
  library(tidyr)
})

# Read PDF filenames from command-line arguments. only use for 2 group comparisons (e.g., Mixed vs Pure; Male vs Female, HSA vs Other)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript /data/rodriguesrr/scripts/R/extract_alpha_div_summary_from_pdfs.R file1.pdf file2.pdf ...")
}
pdf_files <- args

# Define diversity metrics
metrics <- c(
  "Number of Observed Features",
  "Chao1 Index",
  "Shannon Index",
  "Simpson Index",
  "Inverse Simpson Index"
)

# Function to extract from a single PDF
extract_from_pdf <- function(pdf_file) {
  pdf_text_combined <- paste(pdf_text(pdf_file), collapse = "\n")
  
  # Extract group names from legends or capitalized patterns
  legend_matches <- str_match_all(pdf_text_combined, "Legend:?\\s*([A-Za-z0-9_ ]+\\s+vs\\s+[A-Za-z0-9_ ]+)")[[1]]
  if (nrow(legend_matches) > 0) {
    comparison_label <- str_trim(legend_matches[1, 2])
  } else {
    # Fallback to filename
    comparison_label <- tools::file_path_sans_ext(basename(pdf_file))
  }

  # Extract all p-values: handles decimals and scientific notation
  pval_strings <- str_extract_all(pdf_text_combined, "\\b(0|1|0?\\.\\d+|[0-9\\.]+e-?[0-9]+)\\b")[[1]]
  pval_nums <- suppressWarnings(as.numeric(pval_strings))
  pval_nums <- pval_nums[!is.na(pval_nums) & pval_nums >= 0 & pval_nums <= 1]

  if (length(pval_nums) < length(metrics)) {
    warning(paste("Expected", length(metrics), "p-values but found", length(pval_nums), "in", pdf_file, "- skipping."))
    return(NULL)
  }

  # Build data frame
  df <- data.frame(
    Metric = metrics,
    pval_nums[1:length(metrics)]
  )
  colnames(df)[2] <- comparison_label
  return(df)
}

# Run extraction across all files
summary_tables <- lapply(pdf_files, extract_from_pdf)
summary_tables <- summary_tables[!sapply(summary_tables, is.null)]

if (length(summary_tables) == 0) {
  stop("No valid summary tables extracted.")
}

# Merge all by Metric
final_table <- Reduce(function(x, y) full_join(x, y, by = "Metric"), summary_tables)

# Print and save output
print(final_table, row.names = FALSE)

write.csv(final_table, "alpha_diversity_summary.csv", row.names = FALSE)

