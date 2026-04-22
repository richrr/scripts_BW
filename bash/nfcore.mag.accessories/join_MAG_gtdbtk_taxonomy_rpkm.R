#!/usr/bin/env Rscript

# RPKM
  ## log2(RPKM + 1) for heatmap and ordination visualizations
  ## MaAsLin2 (normalization = NONE) for differential abundance testing

# counts
  ## MaAsLin2 (normalization = TSS) for differential abundance testing
  ## DESeq2 or limma-voom for differential abundance testing (use RAW counts, no TSS)

# Method	          Input	    Normalization
# MaAsLin2 (RPKM)	  RPKM	    NONE
# MaAsLin2 (counts)	counts	  TSS (or CLR)
# DESeq2	          counts	  internal
# limma-voom	      counts	  internal

# cd /data/rodriguesrr/Koltsova/analysis/Nov2025_IL22_Alb_Vil/nf-core-mag/odir_mag_no_spades_busco_metabinner
# usage: Rscript /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/join_MAG_gtdbtk_taxonomy_rpkm.R --taxonomy GenomeBinning/DASTool/bins.dRep.MAGs/gtdbtk_r226_out/gtdbtk.bac120.summary.tsv --rpkm GenomeBinning/DASTool/bins.dRep.samples

# cd /data/Trinchieri_lab/rodriguesrr/FMT_nf_mag
# usage: Rscript /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/join_MAG_gtdbtk_taxonomy_rpkm.R --taxonomy DASTool_all/bins.dRep.MAGs/gtdbtk_r226_out/gtdbtk.bac120.summary.tsv,DASTool_all/bins.dRep.MAGs/gtdbtk_r226_out/gtdbtk.ar53.summary.tsv --rpkm DASTool_all/bins.dRep.samples

# --------------------------------------------
# Join GTDB-Tk taxonomy with BBTools RPKM output
# MAG-level abundance via SUM(contig RPKM)
# --------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(readr)
  library(purrr)
  library(stringr)
})

# ---------------------------
# ARGUMENTS
# ---------------------------
option_list <- list(
  make_option("--rpkm", type = "character",
              help = "Directory containing *.rpkm files"),
  make_option("--taxonomy", type = "character",
            help = "Comma-separated GTDB-Tk summary TSVs (bac120,ar53)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$rpkm) || is.null(opt$taxonomy)) {
  stop("ERROR: must provide --rpkm and --taxonomy\n", call. = FALSE)
}

RPKM_DIR <- opt$rpkm

GTDB_FILES <- str_split(opt$taxonomy, ",")[[1]]

print(GTDB_FILES)
# ---- sanity check ----
missing <- GTDB_FILES[!file.exists(GTDB_FILES)]
if (length(missing) > 0) {
  stop("ERROR: taxonomy file(s) not found:\n  ",
       paste(missing, collapse = "\n  "),
       call. = FALSE)
}

gtdb <- map_dfr(GTDB_FILES, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    select(
      mag = user_genome,
      classification
    )
})

head(gtdb)
tail(gtdb)

# ---- check for duplicated MAGs ----
dup <- gtdb$mag[duplicated(gtdb$mag)]
if (length(dup) > 0) {
  stop(
    "ERROR: Duplicate MAGs found across taxonomy files:\n  ",
    paste(unique(dup), collapse = ", "),
    call. = FALSE
  )
}

# ---------------------------
# LOAD + PARSE RPKM FILES
# ---------------------------
EXPECTED_FIELDS <- c(
  "Name", "Length", "Bases", "Coverage",
  "Reads", "RPKM", "Frags", "FPKM"
)

rpkm_files <- list.files(RPKM_DIR, pattern = "\\.rpkm$", full.names = TRUE)

if (length(rpkm_files) == 0) {
  stop("ERROR: no .rpkm files found in ", RPKM_DIR, call. = FALSE)
}

rpkm_long <- map_dfr(rpkm_files, function(f) {

  sample_id <- sub("\\.rpkm$", "", basename(f))

  # ---- Validate header exists (not first line!) ----
  hdr <- read_lines(f, n_max = 25)

header_line <- grep(
  "^\\s*#Name\\tLength\\tBases\\tCoverage\\tReads\\tRPKM\\tFrags\\tFPKM",
  hdr,
  value = TRUE
)

if (length(header_line) == 0) {
  stop(
    paste0(
      "\nERROR: Unexpected .rpkm format in file:\n  ", f,
      "\nExpected header not found:",
      "\n  #Name\\tLength\\tBases\\tCoverage\\tReads\\tRPKM\\tFrags\\tFPKM",
      "\nBBTools output format may have changed."
    ),
    call. = FALSE
  )
}

  # ---- Parse data safely ----
  df <- read_tsv(
    f,
    comment = "#",
    col_names = EXPECTED_FIELDS,
    show_col_types = FALSE
  )

  # ---- Collapse contigs → MAG using SUM ----
  df %>%
    mutate(
      mag = str_remove(Name, "\\$.*$")  # strip contig suffix
    ) %>%
    group_by(mag) %>%
    summarise(
      rpkm = sum(RPKM, na.rm = TRUE),
      counts = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(sample = sample_id)
})

# ---------------------------
# JOIN TAXONOMY
# ---------------------------
mag_abundance_long <- rpkm_long %>%
  left_join(gtdb, by = "mag")


mag_abundance_wide_rpkm <- mag_abundance_long %>%
  select(mag, sample, rpkm) %>%
  pivot_wider(
    names_from  = sample,
    values_from = rpkm,
    values_fill = 0
  ) 


mag_abundance_wide_counts <- mag_abundance_long %>%
  select(mag, sample, counts) %>%
  pivot_wider(
    names_from  = sample,
    values_from = counts,
    values_fill = 0
  ) 


# ---------------------------
# WRITE OUTPUTS
# ---------------------------
write_tsv(
  mag_abundance_long,
  "MAG_sample_rpkm_with_taxonomy.long.tsv"
)

write_tsv(
  mag_abundance_wide_rpkm,
  "MAG_sample_rpkm_with_taxonomy.wide.tsv"
)

write_tsv(
  mag_abundance_wide_counts,
  "MAG_sample_counts_with_taxonomy.wide.tsv"
)

cat("Done.\n")
