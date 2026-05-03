#!/bin/bash
# /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/generate_mag_swarm.sh
# Usage:
# bash /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/generate_mag_swarm.sh <BATCH_DIR> <OUTROOT> -g 400 -t 64 --time 4-00:00:00 --maxrunning 1
# bash /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/generate_mag_swarm.sh batches batches/mag_batches -g 400 -t 64 --time 4-00:00:00 --maxrunning 1
# then: swarm --module nextflow/25.04.2,pandoc --logdir batches/mag_batches/logs -g 400 -t 64 --time 4-00:00:00 --maxrunning 1 ./run-nfcore.mag.swarm
# add --resume tag while rerun.
# or: swarm --module nextflow/25.04.2,pandoc --logdir batches/mag_batches/logs_rerun -g 400 -t 64 --time 4-00:00:00 --maxrunning 1 rerun_failed_batches.swarm

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <BATCH_DIR> <OUTROOT> [-g 240 -t 48 --gres=lscratch:100 --time 4-00:00:00 --maxrunning 1 --sbatch \"--mail-type=ALL --mail-user=rodriguesrr@nih.gov\"]"
    exit 1
fi
RUN_DIR=$(pwd)

BATCH_DIR=$1
OUTROOT=$2
shift 2
SWARM_OPTS="$@"

CONFIG=/data/rodriguesrr/config_files/mag.config
SWARM_FILE="./run-nfcore.mag.swarm"

mkdir -p "$OUTROOT/logs"

# Remove previous swarm file if exists
rm -f "$SWARM_FILE"

# Find all batch CSVs
SAMPLES_BATCHES=("$BATCH_DIR/samples"/*.csv)
ASSEMBLY_BATCHES=("$BATCH_DIR/assembly"/*.csv)

# Make sure we have the same number of sample & assembly files
if [ "${#SAMPLES_BATCHES[@]}" -ne "${#ASSEMBLY_BATCHES[@]}" ]; then
    echo "ERROR: Number of sample CSVs (${#SAMPLES_BATCHES[@]}) != assembly CSVs (${#ASSEMBLY_BATCHES[@]})"
    exit 1
fi

#echo "# Swarm file for nf-core/mag batches" > "$SWARM_FILE"
: > "$SWARM_FILE"

#BATCHES=( $(ls "$BATCH_DIR/samples" | grep -E '^batch[0-9]+\.csv$' | sort) )
mapfile -t BATCHES < <(find "$BATCH_DIR/samples" -maxdepth 1 -name 'batch*.csv' -printf '%f\n' | sort)

# --skip_gtdbtk is done per batch; because we can re-run it once at the end on the non-redunc dRep MAGs
for batch_csv in "${BATCHES[@]}"; do
    BATCH_ID="${batch_csv%.csv}"

    echo "export NXF_WORK='$OUTROOT/nf_work/${BATCH_ID}'; mkdir -p \$NXF_WORK; \
export NXF_JVM_ARGS='-Xms8g -Xmx64g'; \
export NXF_SINGULARITY_CACHEDIR=/data/\$USER/nxf_singularity_cache; \
export SINGULARITY_CACHEDIR=/data/\$USER/.singularity; \
nextflow run nf-core/mag -r 5.3.0 -profile biowulf \
--input '$BATCH_DIR/samples/${BATCH_ID}.csv' \
--assembly_input '$BATCH_DIR/assembly/${BATCH_ID}.csv' \
--outdir '$OUTROOT/${BATCH_ID}' \
--skip_gtdbtk --gtdb_db /data/rodriguesrr/db/gtdbtk_r226/gtdbtk_r226_data.tar.gz \
--refine_bins_dastool \
--binning_map_mode group \
--postbinning_input refined_bins_only \
--run_checkm2 \
--checkm2_db /data/rodriguesrr/db/checkm2.v3/uniref100.KO.1.dmnd \
--skip_megahit \
--skip_spades \
--skip_prodigal \
--skip_prokka \
--skip_metaeuk \
--skip_concoct \
--skip_comebin \
--skip_metabinner \
--skip_semibin \
--run_busco false \
--run_checkm false \
-c $CONFIG" >> "$SWARM_FILE"

done


echo "Swarm file generated: $SWARM_FILE"
echo "Note that the sbatch args must be in the double quotes when actually launching the swarm."
echo "Launch with: swarm --module nextflow/25.04.2,pandoc --logdir $OUTROOT/logs $SWARM_OPTS $SWARM_FILE"

