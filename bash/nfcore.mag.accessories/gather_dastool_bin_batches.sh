#!/bin/bash
set -euo pipefail

# cd /data/Trinchieri_lab/rodriguesrr/FMT_nf_mag
# usage: bash /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/gather_dastool_bin_batches.sh batches/mag_batches /data/Trinchieri_lab/rodriguesrr/FMT_nf_mag/DASTool_all

# using the checkM2 QC and only taking good fasta files from DASTool bins from batches in $MAG_ROOT
# does not use raw dastool bins, unbinned or low quality MAGs.

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <MAG_ROOT> <FINAL_OUT>"
    echo "Example: $0 /path/to/nfcore_mag_batches /path/to/combined_dastool"
    exit 1
fi

MAG_ROOT=$1
FINAL_OUT=$2

# ---------------------------
# QC thresholds (standard MAG cutoffs)
# ---------------------------
MIN_COMPLETENESS=50
MAX_CONTAMINATION=10


mkdir -p "$FINAL_OUT/bins"

COMBINED_QC="$FINAL_OUT/checkm2_filtered.tsv"
GENOME_LIST="$FINAL_OUT/genomes_list.txt"

echo "Gathering script is using the checkM2 QC and only taking good fasta files from DASTool bins from batches in $MAG_ROOT ..."
#echo "Final output will go to $FINAL_OUT"
echo "MAG root:   $MAG_ROOT"
echo "Output:     $FINAL_OUT"
echo "QC: comp >= $MIN_COMPLETENESS | cont <= $MAX_CONTAMINATION"


# ---------------------------
# Temporary associative arrays (bash 4+)
# ---------------------------
declare -A QC_COMP
declare -A QC_CONT

echo -e "Batch\tName\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tTotal_Contigs\tMax_Contig_Length\tAdditional_Notes" > "$COMBINED_QC"


# ---------------------------
# STEP 1: Parse ALL CheckM2 results
# the following code parses the checkM2 summary files from each batch, applies the QC thresholds, and builds a combined QC table with batch info. 
# This allows us to keep track of which bins passed QC and their associated metadata for downstream filtering and analysis.
# ---------------------------
echo "Parsing CheckM2 files..."

for BATCH_DIR in "$MAG_ROOT"/batch*; do
    BATCH=$(basename "$BATCH_DIR")
    QC_FILE="$BATCH_DIR/GenomeBinning/QC/checkm2_summary.tsv"

    if [[ -f "$QC_FILE" ]]; then
        tail -n +2 "$QC_FILE" | awk -v b="$BATCH" -v minc="$MIN_COMPLETENESS" -v maxc="$MAX_CONTAMINATION" '
        BEGIN {FS=OFS="\t"}
        {
            name=$1
            comp=$2
            cont=$3

            if (comp >= minc && cont <= maxc) {
                print b "\t" $0
            }
        }' >> "$COMBINED_QC"
    else
        echo "WARNING: missing $QC_FILE"
    fi
done

# ---------------------------
# STEP 2: Build allowed genome set
# ---------------------------
echo "Building allowed genome list..."

while IFS=$'\t' read -r col1 NAME rest; do
    QC_COMP["$NAME"]=1
done < <(tail -n +2 "$COMBINED_QC")

echo "QC entries loaded: ${#QC_COMP[@]}"


# echo "Lets see the first 10 QC-passing bins:"
# # echo the first 10 QC-passing bins without using head or printf since it crashes with set -e
# count=0
# for name in "${!QC_COMP[@]}"; do
#     echo "$name"
#     count=$((count + 1))
#     if [[ $count -ge 10 ]]; then
#         break
#     fi
# done


# ---------------------------
# STEP 3: Collect bins (ONLY QC-passing)
# ---------------------------
echo "Collecting bins..."

for BATCH_DIR in "$MAG_ROOT"/batch*; do
    BIN_DIR="$BATCH_DIR/GenomeBinning/DASTool/bins"

    if [[ -d "$BIN_DIR" ]]; then
        for f in "$BIN_DIR"/*.fa; do
            [[ -e "$f" ]] || continue
            #echo "Processing $f"
            base=$(basename "$f" .fa)
            #echo "Base name: $base"
            # only keep if QC-passing
            if [[ -n "${QC_COMP[$base]+x}" ]]; then
                ln -sf "$(realpath "$f")" "$FINAL_OUT/bins/"
                #echo "Added $base to final bins"
            fi
        done
    fi
done


# ---------------------------
# STEP 4: Generate dRep genome list
# ---------------------------
echo "Generating dRep genome list..."

# how many bins exist in folder (including symlinks)
find "$FINAL_OUT/bins" -maxdepth 1 -name "*.fa" > "$GENOME_LIST"

# how many symlinks point to .fa files
#find -L "$FINAL_OUT/bins" -maxdepth 1 -type f -name "*.fa" > "$GENOME_LIST"

NUM=$(wc -l < "$GENOME_LIST")

echo "QC-passing genomes: $NUM"

if [[ "$NUM" -eq 0 ]]; then
    echo "ERROR: No QC-passing genomes found"
    exit 1
fi


echo "Bins:        $FINAL_OUT/bins"
echo "QC table:    $COMBINED_QC"
echo "Genome list: $GENOME_LIST"

echo "Done."