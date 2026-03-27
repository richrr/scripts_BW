#!/bin/bash
set -euo pipefail

# cd /data/Trinchieri_lab/rodriguesrr/FMT_nf_mag
# usage: bash /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/gather_dastool_bin_batches.sh batches/mag_batches /data/Trinchieri_lab/rodriguesrr/FMT_nf_mag/DASTool_all


if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <MAG_ROOT> <FINAL_OUT>"
    echo "Example: $0 /path/to/nfcore_mag_batches /path/to/combined_dastool"
    exit 1
fi

MAG_ROOT=$1
FINAL_OUT=$2

mkdir -p "$FINAL_OUT/bins"
#mkdir -p "$FINAL_OUT/logs"

echo "Gathering DASTool bins from batches in $MAG_ROOT ..."
echo "Final output will go to $FINAL_OUT"

# Loop over all batch directories
for BATCH_DIR in "$MAG_ROOT"/batch*/GenomeBinning/DASTool/bins; do
    if [[ -d "$BATCH_DIR" ]]; then
        echo "Processing $BATCH_DIR ..."
        #cp -r "$BATCH_DIR"/* "$FINAL_OUT/bins/"
        for f in "$BATCH_DIR"/*; do
            ln -s "$(realpath "$f")" "$FINAL_OUT/bins/"
        done
    else
        echo "WARNING: No DASTool bins found in $BATCH_DIR"
    fi
done

# Optional: gather logs if they exist
#for LOG_DIR in "$MAG_ROOT"/batch*/logs; do
#    if [[ -d "$LOG_DIR" ]]; then
#        cp -r "$LOG_DIR"/* "$FINAL_OUT/logs/"
#    fi
#done

echo "All DASTool bins (and maybe logs) have been symlinked to $FINAL_OUT"
