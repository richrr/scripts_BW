#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tools)
})

#ml R/4.5.0
#ml python/3.10 

#cd /data/MicrobiomeCore/analysis/ADAD02_20250702/combined_ADAD_01_02_T1/T1_survival_analysis
# Rscript /data/rodriguesrr/scripts/R/prepare_inputs_for_survival_analysis.R patient_data ../../map-files/combined_ADAD_01_02_map_T1/final_T1_sequencing_ids.tsv SampleID "Survival Time (days)" "Dead (1) or Alive/Unknown (0)"

#cd /data/MicrobiomeCore/analysis/ADAD01_20240215/survival_analysis
#Rscript /data/rodriguesrr/scripts/R/prepare_inputs_for_survival_analysis.R test_wrapper_code ../map-files/T1-metadata-for-use_w_Abx_info_Sep2025.txt SampleID "Surv_days" "survival_0_alive_unknown__1_dead"


# ---------------------------
# Parse arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript prepare_survival_inputs.R <prefix> <infile> <sampleID_col_or_absent> <timecol[,unit]> <eventcol> [min_samples]")
}

prefix          <- args[1]
infile          <- args[2]
sample_col_arg  <- args[3]
time_info       <- args[4]
event_col       <- args[5]
min_samples     <- ifelse(length(args) >= 6, args[6], NA)


# ---------------------------
# Read input
# ---------------------------
data <- fread(infile)
cat("Loaded", nrow(data), "rows and", ncol(data), "columns from", infile, "\n")
#head(data)


# ---------------------------
# Extract survival columns
# ---------------------------
time_split <- strsplit(time_info, ",")[[1]]

# Handle cases like "Survival Time (days)" or "OS (months)"
if (length(time_split) == 2) {
  time_col <- time_split[1]
  time_unit <- tolower(time_split[2])
  if (!time_col %in% names(data)) {
  stop(paste("Time column", time_col, "not found in data"))
  }
} else {
  time_col <- time_info
  if (!time_col %in% names(data)) {
  stop(paste("Time column", time_col, "not found in data"))
  }

  # Try to detect unit from the column name itself
  if (grepl("day", time_col, ignore.case = TRUE)) {
    time_unit <- "days"
  } else if (grepl("month", time_col, ignore.case = TRUE)) {
    time_unit <- "months"
  } else if (grepl("year", time_col, ignore.case = TRUE)) {
    time_unit <- "years"
  } else {
    time_unit <- "months" # default
    warning("No explicit time unit found — assuming months.")
  }
}


if (!event_col %in% names(data)) {
  stop(paste("Event column", event_col, "not found in data"))
}

# ---------------------------
# Convert survival time to months
# ---------------------------
convert_to_months <- function(x, unit) {
  if (unit %in% c("day", "days")) x / 30
  else if (unit %in% c("year", "years")) x * 12
  else if (unit %in% c("month", "months")) x
  else stop("Unknown time unit. Use days, months, or years.")
}

data$dmfs_time <- convert_to_months(data[[time_col]], time_unit)
data$dmfs_event <- as.numeric(data[[event_col]])


# ---------------------------
# Remove sample ID column if present
# ---------------------------
if (tolower(sample_col_arg) != "absent" && sample_col_arg %in% names(data)) {
  data <- data %>% select(-all_of(sample_col_arg))
}

# Remove original time/event columns
data <- data %>% select(-all_of(c(time_col, event_col)))
#head(data)


# ---------------------------
# Detect column types
# ---------------------------
is_date <- function(x) {
  if (inherits(x, "Date")) return(TRUE)
  x_chr <- na.omit(as.character(x))
  if (length(x_chr) == 0) return(FALSE)
  all(grepl("^\\d{4}-\\d{2}-\\d{2}$", x_chr) | 
      grepl("^\\d{1,2}/\\d{1,2}/\\d{2,4}$", x_chr))
}

is_binary <- function(x) {
  x_num <- suppressWarnings(as.numeric(na.omit(x)))
  ux <- unique(x_num)
  all(ux %in% c(0, 1)) && length(ux) <= 2
}

categorical_cols <- c()
continuous_cols <- c()
removed_cols <- c()

for (colname in setdiff(names(data), c("dmfs_time", "dmfs_event"))) {
  x <- data[[colname]]
  
  if (is_date(x)) {
    next
  } else if (is_binary(x)) {
    x_num <- suppressWarnings(as.numeric(x))
    data[[colname]] <- ifelse(is.na(x_num), NA,
                              ifelse(x_num == 1, "yes", "no"))
    categorical_cols <- c(categorical_cols, colname)
  } else if (is.character(x) || is.factor(x)) {
    n_levels <- length(unique(na.omit(x)))
    
    if (n_levels > 5) {
      removed_cols <- c(removed_cols, colname)
      next  # skip this column entirely
    } else {
      categorical_cols <- c(categorical_cols, colname)
      if (n_levels > 3) {
        warning(sprintf("Column '%s' has %d levels (more than 3)", colname, n_levels))
      }
    }
  } else if (is.numeric(x)) {
    continuous_cols <- c(continuous_cols, colname)
  }
}

if (length(removed_cols) > 0) {
  cat("\n\nRemoved categorical columns with >5 levels:\n  ",
      paste(removed_cols, collapse = ", "), "\n\n\n")
}


cat("\n\nCategorical columns:", paste(categorical_cols, collapse=", "), "\n")
cat("\n\nContinuous columns:", paste(continuous_cols, collapse=", "), "\n")

# ---------------------------
# Prepare output directories
# ---------------------------
dir_cont <- paste0(prefix, "_continuous")
dir_cat  <- paste0(prefix, "_categories")
dir.create(dir_cont, showWarnings = FALSE)
dir.create(dir_cat, showWarnings = FALSE)

# ---------------------------
# Output files
# ---------------------------
out_cont_file <- file.path(dir_cont, paste0(prefix, "_continuous.txt"))
out_cat_file  <- file.path(dir_cat,  paste0(prefix, "_categories.txt"))

# For continuous data
cont_df <- data %>%
  select(all_of(c(continuous_cols, "dmfs_time", "dmfs_event"))) %>%
  distinct()

fwrite(cont_df, out_cont_file, sep="\t", quote=FALSE, na="NA")
cat("\nSaved continuous file to:", out_cont_file, "\n")

# For categorical data
cat_df <- data %>%
  select(all_of(c(categorical_cols, "dmfs_time", "dmfs_event"))) %>%
  distinct()

fwrite(cat_df, out_cat_file, sep="\t", quote=FALSE, na="NA")
cat("Saved categorical file to:", out_cat_file, "\n")

# ---------------------------
# Print next-step messages
# ---------------------------
cat("\nNext steps:\n")
cat("\nContinuous data:\n")
if (!is.na(min_samples)) {
  cat("  Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/EvaluateCutpoints_multiple-analyses_analysis.R",
      dir_cont, paste0(prefix, "_continuous.txt"), min_samples, "\n")
} else {
  cat("  Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/EvaluateCutpoints_multiple-analyses_analysis.R",
      dir_cont, paste0(prefix, "_continuous.txt"), "\n")
}
cat("\nCategorical data:\n")
cat("  Rscript /data/rodriguesrr/scripts/R/Evaluate_cutpoints/Cox_Regression_Plot_Kaplan_Meyer_for_categories_multiple-analyses.R",
    dir_cat, paste0(prefix, "_categories.txt") , "\n")

