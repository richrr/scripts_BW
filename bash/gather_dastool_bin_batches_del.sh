#!/bin/bash
# gather_dastool_bin_batches.sh
# Usage: bash /data/rodriguesrr/scripts/bash/gather_dastool_bin_batches.sh <DASTOOL_PARENT_DIR> <OUTPUT_DIR>

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <DASTOOL_PARENT_DIR> <OUTPUT_DIR>"
    exit 1
fi

DASTOOL_PARENT="$1"
OUTDIR="$2"

mkdir -p "$OUTDIR"

echo "Combining DASTool-refined bins from batches in: $DASTOOL_PARENT"
echo "Output will go to: $OUTDIR"

# Loop through batch directories
for batch_dir in "$DASTOOL_PARENT"/*/; do
    batch_name=$(basename "$batch_dir")
    echo "Processing batch: $batch_name"

    # Find all refined bins
    for bin in "$batch_dir"/*.fa; do
        bin_name=$(basename "$bin")
        # Prefix with batch to avoid collisions
        cp "$bin" "$OUTDIR/${batch_name}_${bin_name}"
    done
done

echo "All bins combined into: $OUTDIR"
echo "You can now run dRep on $OUTDIR/*.fa"

