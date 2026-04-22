### bbtools_indexing_read_quantification ###
# Input
#     dereplicated MAGs
#     sample sheet pointing to cleaned reads
# Output
#     MAG × sample abundance table

# cd /data/Trinchieri_lab/rodriguesrr/FMT_nf_mag
# usage: bash /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/bbtools_indexing_read_quantification.sh /data/Trinchieri_lab/rodriguesrr/FMT_nf_mag/DASTool_all/bins fmt-samples-w-jamsCleanReads.csv
# swarm --logdir logs -g 300 -t 32 --time 8-00:00:00 run_bbtools_index_quant.sh

#cd /data/rodriguesrr/Koltsova/analysis/Nov2025_IL22_Alb_Vil/nf-core-mag/odir_mag_no_spades_busco_metabinner
#bash /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/bbtools_indexing_read_quantification.sh GenomeBinning/DASTool/bins ../samples-w-CleanReads.csv
#swarm --logdir swarmlogs -g 300 -t 32 --time 2-00:00:00 run_bbtools_index_quant.sh


#!/usr/bin/env bash

set -euo pipefail

# ---------------------------
# USAGE
# ---------------------------
if [ $# -ne 2 ]; then
    echo "Usage: $0 <bins_directory> <Three_column_sample_sheet.csv>"
    exit 1
fi

BINS_DIR=$1            # e.g., DASTool_all/bins # DO NOT GIVE IT DASTool_all/bins.dRep.MAGs OR DASTool_all/bins.dRep.MAGs/dereplicated_genomes/
SAMPLE_SHEET=$2        # file with sample names, and paths to cleaned R1 and R2


# ---------------------------
# NORMALIZE SAMPLE SHEET
# ---------------------------
# Ensure last line ends with a newline and remove Windows CRs if present

# Check last character
last_char=$(tail -c1 "$SAMPLE_SHEET" || true)

if [[ "$last_char" != $'\n' && "$last_char" != $'\r' ]]; then
    echo "Normalizing sample sheet: adding missing trailing newline"
    sed -i '$a\' "$SAMPLE_SHEET"
fi

# Convert any CRLF (Windows) endings to LF (Unix)
sed -i 's/\r$//' "$SAMPLE_SHEET"



# ---------------------------
# OUTPUT DIRECTORIES
# ---------------------------
DREP_OUT="${BINS_DIR}.dRep.MAGs"
SAMPLE_OUT="${BINS_DIR}.dRep.samples"

mkdir -p "$DREP_OUT" "$SAMPLE_OUT"

echo "Bins directory:           $BINS_DIR"
echo "Dereplicated bins out:    $DREP_OUT"
echo "Sample output directory:  $SAMPLE_OUT"
echo "Sample sheet:             $SAMPLE_SHEET"

# ---------------------------
# LOAD MODULES
# ---------------------------

module load bbtools/39.06


# ---------------------------
# INDEXING
# ---------------------------
REF_INDEX="$DREP_OUT/dereplicated_genomes.ref"

if [[ -d "$REF_INDEX" ]]; then
    echo "BBTools index already exists → skipping indexing."
else
    echo "Indexing dereplicated MAGs..."
    unset _JAVA_OPTIONS # prevent Java heap conflicts
    bbtools bbsplit -Xmx64g \
      ref="$DREP_OUT/dereplicated_genomes/" \
      path="$REF_INDEX"
fi


# ---------------------------
# BUILD SWARM FILE (PER-SAMPLE BBTOOLS)
# ---------------------------
echo "Building swarm file for per-sample quantification..."

SWARM_FILE="$SAMPLE_OUT/bbtools_quant.swarm"
rm -f "$SWARM_FILE"

TOTAL=0
SKIPPED_DONE=0
SKIPPED_MISSING=0
ADDED=0

# Read CSV, skip header, extract first three columns
#tail -n +2 "$SAMPLE_SHEET" | while IFS=',' read -r SAMPLE R1 R2 || [ -n "$SAMPLE" ]; do
TMP_SHEET="$SAMPLE_OUT/tmp_samples.csv"
tail -n +2 "$SAMPLE_SHEET" > "$TMP_SHEET"

while IFS=',' read -r SAMPLE R1 R2 || [ -n "$SAMPLE" ]; do
    # Trim whitespace
    SAMPLE=$(echo "$SAMPLE" | tr -d '[:space:]')
    R1=$(echo "$R1" | tr -d '[:space:]')
    R2=$(echo "$R2" | tr -d '[:space:]')

    if [[ -z "$SAMPLE" ]]; then
        continue
    fi

    TOTAL=$((TOTAL + 1))

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "WARNING: FASTQs not found for $SAMPLE  → skipping."
        echo "  R1: $R1"
        echo "  R2: $R2"
        SKIPPED_MISSING=$((SKIPPED_MISSING + 1))
        continue
    fi

    RPKM_OUT="$SAMPLE_OUT/${SAMPLE}.rpkm"
    LOG="$SAMPLE_OUT/${SAMPLE}.bbsplit.log"

    # Skip already processed samples (resume-safe)
    if [[ -s "$RPKM_OUT" ]]; then
        echo "Sample $SAMPLE already quantified → skipping."
        SKIPPED_DONE=$((SKIPPED_DONE + 1))
        continue
    fi

    
    # echo "Processing sample: $SAMPLE"
    # echo "  R1: $R1"
    # echo "  R2: $R2"

    # unset _JAVA_OPTIONS
    # LOG="$SAMPLE_OUT/${SAMPLE}.bbsplit.log"

    # bbtools bbsplit -Xmx64g \
    #   unpigz=t \
    #   threads=8 \
    #   minid=0.90 \
    #   sortscafs=f \
    #   nzo=f \
    #   path="$REF_INDEX" \
    #   statsfile="$SAMPLE_OUT/${SAMPLE}.statsfile" \
    #   scafstats="$SAMPLE_OUT/${SAMPLE}.scafstats" \
    #   covstats="$SAMPLE_OUT/${SAMPLE}.covstat" \
    #   rpkm="$SAMPLE_OUT/${SAMPLE}.rpkm" \
    #   refstats="$SAMPLE_OUT/${SAMPLE}.refstats" \
    #   in="$R1" \
    #   in2="$R2" &> "$LOG"

    echo "Adding $SAMPLE to swarm..."

    echo "unset _JAVA_OPTIONS; \
    bbtools bbsplit -Xmx64g \
        unpigz=t \
        threads=16 \
        minid=0.90 \
        sortscafs=f \
        nzo=f \
        path=\"$REF_INDEX\" \
        statsfile=\"$SAMPLE_OUT/${SAMPLE}.statsfile\" \
        scafstats=\"$SAMPLE_OUT/${SAMPLE}.scafstats\" \
        covstats=\"$SAMPLE_OUT/${SAMPLE}.covstat\" \
        rpkm=\"$RPKM_OUT\" \
        refstats=\"$SAMPLE_OUT/${SAMPLE}.refstats\" \
        in=\"$R1\" \
        in2=\"$R2\" &> \"$LOG\"" >> "$SWARM_FILE"
    ADDED=$((ADDED + 1))
done < "$TMP_SHEET"

rm -f "$TMP_SHEET"


# # Verify per-sample outputs at the very end
# #EXPECTED=$(tail -n +2 "$SAMPLE_SHEET" | wc -l)
# EXPECTED=$(tail -n +2 "$SAMPLE_SHEET" | awk 'NF {count++} END {print count}')
# OBSERVED=$(ls -1 "$SAMPLE_OUT"/*.rpkm 2>/dev/null | wc -l)
# echo "Expected samples: $EXPECTED"
# echo "Observed RPKM files: $OBSERVED"

# if [[ "$OBSERVED" -ne "$EXPECTED" ]]; then
#     echo "WARNING: Not all samples were successfully quantified."
#     echo "You can safely rerun the job to resume."
# else
#     echo "All samples successfully quantified."
# fi

#NUM_JOBS=$(wc -l < "$SWARM_FILE")

# ---------------------------
# SUMMARY
# ---------------------------
echo "----------------------------------------"
echo "Swarm file: $SWARM_FILE"
echo ""
echo "Total samples in sheet:   $TOTAL"
echo "Jobs added to swarm:      $ADDED"
echo "Skipped (already done):   $SKIPPED_DONE"
echo "Skipped (missing FASTQ):  $SKIPPED_MISSING"
echo ""
echo "Run with:"
echo "swarm --module=bbtools/39.06 -f $SWARM_FILE -g 150 -t 16 --time 1-00:00:00 --maxrunning 15 --logdir logs"
echo "----------------------------------------"

echo "NOTE: Verification block should be run AFTER swarm completes."

echo "DONE."
