#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
})

# Path to your folder
#indir <- "path/to/your/files"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript summarize_clusterprofiler_results.R <indir>")
}
indir <- args[1]


# Databases list
db_list <- c("Cell_marker_2.0_Human", "REACTOME", "GOBP", "KEGG_LEGACY", 
             "KEGG_MEDICUS", "HALLMARK", "WIKIPATHWAYS", "BIOCARTA", 
             "IMMUNESIGDB", "CellTypeSig")

# Escape regex characters in db names for pattern matching
db_pattern <- paste0("\\b(", paste0(str_replace_all(db_list, "\\.", "\\\\."), collapse = "|"), ")\\b")

# List all csv files
files <- list.files(indir, pattern = "\\.csv$", full.names = TRUE)

read_with_meta <- function(f) {
  fname <- basename(f)
  
  # Extract comparison-condition (before first .csv_)
  comparison_condition <- str_extract(fname, "^[^_]+_?[^_]+(?=\\.csv)")
  
  # Extract database (between .csv_ and _<direction>_ORA.csv)
  database <- str_match(fname, "\\.csv_(.*)_(up|down|all)_ORA\\.csv$")[,2]
  
  # Extract direction (up|down|all)
  direction <- str_match(fname, "_(up|down|all)_ORA\\.csv$")[,2]
  
  # Split comparison-condition into two parts
  parts <- str_split(comparison_condition, "_", n = 2)[[1]]
  comparison <- parts[1]
  condition  <- parts[2]
  
  
  # Read CSV (don't auto-guess types too hard)
  dat <- suppressMessages(read_csv(f, col_types = cols(.default = "c")))
  
  # Coerce numeric columns
  numeric_cols <- c("RichFactor", "FoldEnrichment", "zScore", 
                    "pvalue", "p.adjust", "qvalue", "Count")
  dat <- dat %>%
    mutate(across(any_of(numeric_cols), ~ suppressWarnings(as.numeric(.x))))
  
  # Add metadata
  dat %>%
    mutate(
      comparison = comparison,
      condition  = condition,
      database   = database,
      direction  = direction,
      sourcefile = fname
    )
}

# Read and combine
all_data <- map_dfr(files, read_with_meta)

# Save
write_csv(all_data, file.path(indir, "All_ORA_results_combined.csv"))


