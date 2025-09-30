args = commandArgs(trailingOnly=TRUE)

# usage: ml R/4.4.3
# cd /data/MGC_CCRSF/testmedi_pittsdata/JAMS_plots/Transnet
# Rscript /data/rodriguesrr/scripts/R/parse-and-summarize-comp-results.R analysis/food__content/mw_sp_comp-output.csv analysis/food_subgroup__abundance/mw_sp_comp-output.csv analysis/food_group__abundance/mw_sp_comp-output.csv  analysis/wikipedia_id__abundance/mw_sp_comp-output.csv

for(infile in args)
{
    print(infile)

    df = read.csv(infile, header=T, row.names=1, check.names=F)

    # Identify columns
    fc_cols <- grep("FolChMedian", names(df), value = TRUE)
    pval_cols <- grep("pvalue", names(df), value = TRUE)

    # Convert p-value columns to numeric (if needed)
    df[pval_cols] <- lapply(df[pval_cols], function(x) as.numeric(as.character(x)))

    # Select rows with any p-value < 0.1
    rows_with_low_pval <- apply(df[pval_cols], 1, function(row) any(row < 0.1, na.rm = TRUE))

    #print(length(rows_with_low_pval))

    # Subset the data frame
    paired_cols <- grep("FolChMedian|pvalue", colnames(df), value = TRUE)
    df_filtered <- df[rows_with_low_pval, paired_cols]

    # Split into fold change and p-value columns (alternating)
    fc_cols <- paired_cols[seq(1, length(paired_cols), by = 2)]
    pval_cols <- paired_cols[seq(2, length(paired_cols), by = 2)]

    # Initialize counters
    pos_fc_count <- numeric(nrow(df_filtered))
    neg_fc_count <- numeric(nrow(df_filtered))

    # Count positive and negative fold changes per row if its pvalue is signif
    for (i in seq_along(fc_cols)) {
    fc_vals   <- as.numeric(df_filtered[[fc_cols[i]]])
    pval_vals <- as.numeric(df_filtered[[pval_cols[i]]])
    
    pos_fc_count <- pos_fc_count + ifelse(!is.na(pval_vals) & pval_vals < 0.1 & fc_vals > 0, 1, 0)
    neg_fc_count <- neg_fc_count + ifelse(!is.na(pval_vals) & pval_vals < 0.1 & fc_vals < 0, 1, 0)
    }


    # Add those as new columns
    df_filtered$pos_fc_count <- pos_fc_count
    df_filtered$neg_fc_count <- neg_fc_count


    write.csv(df_filtered, paste0(infile, ".summary.csv"), row.names=T)
}