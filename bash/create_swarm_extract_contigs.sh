#!/usr/bin/env bash
set -euo pipefail

# cd /data/Trinchieri_lab/rodriguesrr/cohorts_nf_mag_NatMed2022/JAMStarballs
# usage: bash /data/rodriguesrr/scripts/bash/create_swarm_extract_contigs.sh ../arranged_155_samples.txt

# Defaults
DEFAULT_CSV="../sample-list.txt"
DEFAULT_OUTDIR="extracted"
DEFAULT_SWARM="extract_contigs.swarm"

# User inputs (fallback to defaults if not provided)
CSV="${1:-$DEFAULT_CSV}"
OUTDIR="${2:-$DEFAULT_OUTDIR}"
SWARM="${3:-$DEFAULT_SWARM}"

echo "Using:"
echo "  CSV    = $CSV"
echo "  OUTDIR = $OUTDIR"
echo "  SWARM  = $SWARM"


# Step 1: create output directory
mkdir -p "$OUTDIR"
: > "$SWARM"

# Step 2: read only the first column (sample), skip header
tail -n +2 "$CSV" | cut -d',' -f1 | while IFS= read -r sample; do
  tgz="${sample}_JAMS.tar.gz"
  out_fa="${OUTDIR}/${sample}_JAMS/${sample}_contigs.fasta"

  # Skip if tar.gz does not exist
  [[ -f "$tgz" ]] || continue

  # skip if contigs already extracted
  [[ -f "$out_fa" ]] && continue

  cat >> "$SWARM" <<EOF
tar -xzf "$tgz" -C "$OUTDIR" "${sample}_JAMS/${sample}_contigs.fasta"
EOF
done

echo "Wrote swarm file: $SWARM"
