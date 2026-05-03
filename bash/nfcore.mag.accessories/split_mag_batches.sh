#!/bin/bash
set -euo pipefail

# --** make sure the group column for all samples is the same in the samplesheet and assembly input files. 
# also make sure there is an empty line at the end. the last sample line should have an next line character after it.**--
# the script splits these files into smaller batches. 
# The script will then overwrite that column for each batch when writing the per-batch CSVs (so batch01, batch02, etc.).
# Usage: bash /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/split_mag_batches.sh fmt-samples-w-jamsContigs.csv assembly_input-w-jamsContigs.csv batches 35

SAMPLESHEET=$1
ASSEMBLY_SHEET=$2
OUTDIR=$3
BATCH_SIZE=${4:-35} # “Set BATCH_SIZE to the 4th argument of the script, but if no 4th argument is given, default to 35 samples per batch.”

mkdir -p "$OUTDIR"/{samples,assembly}

# Extract headers
head -n 1 "$SAMPLESHEET" > "$OUTDIR/header.samples.csv"
head -n 1 "$ASSEMBLY_SHEET" > "$OUTDIR/header.assembly.csv"

# Remove headers, keep body
tail -n +2 "$SAMPLESHEET" > "$OUTDIR/samples.body.csv"
tail -n +2 "$ASSEMBLY_SHEET" > "$OUTDIR/assembly.body.csv"

TOTAL=$(wc -l < "$OUTDIR/samples.body.csv")
NBATCHES=$(( (TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))

echo "Total samples: $TOTAL"
echo "Batch size:    $BATCH_SIZE"
echo "Batches:       $NBATCHES"

# For each batch:
    # Pick the rows for this batch (sed).
    # Set the group column to the batch ID (awk).
    # Prepend the header (cat).
    # Save to a new CSV (>).

# Split into batches
for ((i=0; i<NBATCHES; i++)); do
  BATCH_ID=$(printf "batch%02d" $((i+1)))
  START=$(( i * BATCH_SIZE + 1 ))
  END=$(( (i+1) * BATCH_SIZE ))
  
  # Samplesheet
  sed -n "${START},${END}p" "$OUTDIR/samples.body.csv" \
    | awk -F',' -v OFS=',' -v b="$BATCH_ID" '{$2=b; print}' \
    | cat "$OUTDIR/header.samples.csv" - \
    > "$OUTDIR/samples/${BATCH_ID}.csv"

  # Assembly input
  sed -n "${START},${END}p" "$OUTDIR/assembly.body.csv" \
    | awk -F',' -v OFS=',' -v b="$BATCH_ID" '{$2=b; print}' \
    | cat "$OUTDIR/header.assembly.csv" - \
    > "$OUTDIR/assembly/${BATCH_ID}.csv"

done

echo "Batching complete."
echo


echo "Run the batches with Swarm:"
echo "  bash /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/generate_mag_swarm.sh $OUTDIR $OUTDIR/mag_batches -g 240 -t 48 --gres=lscratch:100 --time 4-00:00:00 --maxrunning 1 --sbatch \"--mail-type=ALL --mail-user=rodriguesrr@nih.gov\""
echo
echo "Then launch:"
echo "  swarm --module nextflow/25.04.2,pandoc --logdir $OUTDIR/mag_batches/logs -g 240 -t 48 --gres=lscratch:100 --time 4-00:00:00 --maxrunning 1 --sbatch \"--mail-type=ALL --mail-user=rodriguesrr@nih.gov\" ./run-nfcore.mag.swarm "
