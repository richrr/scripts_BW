### dRep and gtdbtk ###
# Input (No reads involved)
#   bins.dRep.MAGs/dereplicated_genomes/*.fa   ← FINAL MAGs
# Output
#   Run GTDB-Tk on dereplicated_genomes/ & produces 
#   bins.dRep.MAGs/gtdbtk_r226_out/
#   ├── classify/
#   │   ├── gtdbtk.bac120.summary.tsv
#   │   ├── gtdbtk.bac120.classification_pplacer.tsv
#   │   └── ...
#   ├── logs/
#   └── gtdbtk.log


# cd /data/Trinchieri_lab/rodriguesrr/FMT_nf_mag
# usage: bash /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/dRep_MAGs_gtdbtk_taxonomy.sh /data/Trinchieri_lab/rodriguesrr/FMT_nf_mag/DASTool_all/bins
# swarm --logdir logs -g 300 -t 32 --time 8-00:00:00 run_dRep_gtdbtk.sh

#cd /data/rodriguesrr/Koltsova/analysis/Nov2025_IL22_Alb_Vil/nf-core-mag/odir_mag_no_spades_busco_metabinner
#bash /data/rodriguesrr/scripts/bash/nfcore.mag.accessories/dRep_MAGs_gtdbtk_taxonomy.sh GenomeBinning/DASTool/bins
#swarm --logdir swarmlogs -g 300 -t 32 --time 2-00:00:00 run_dRep_gtdbtk.sh

#!/usr/bin/env bash

set -euo pipefail

# ---------------------------
# USAGE
# ---------------------------
if [ $# -ne 1 ]; then
    echo "Usage: $0 <bins_directory>"
    exit 1
fi

BINS_DIR=$1            # e.g., DASTool_all/bins

# ---------------------------
# OUTPUT DIRECTORIES
# ---------------------------
DREP_OUT="${BINS_DIR}.dRep.MAGs"

mkdir -p "$DREP_OUT"

echo "Bins directory:           $BINS_DIR"
echo "Dereplicated bins out:    $DREP_OUT"

# ---------------------------
# LOAD MODULES
# ---------------------------
# tell Lmod not to fail the script if a module complains
#set +e
module load drep/3.4.2
#set -e


#module load gtdb-tk/2.3.2
module load gtdb-tk/2.6.1 # uses ref db v226 by default


# ---------------------------
# DEREPLICATION
# ---------------------------
echo "Counting bins..."
#NUM_BINS=$(ls -1 "$BINS_DIR"/*.fa 2>/dev/null | wc -l) # failed due to 126:0 error which is possibly 
                                                                                # because glob+ls takes few seconds if there are 60K fa files
GENOME_LIST="$DREP_OUT/genomes_list.txt"
find -L "$BINS_DIR" -maxdepth 1 -type f -name "*.fa" > "$GENOME_LIST"
NUM_BINS=$(wc -l < "$GENOME_LIST")

if [ "$NUM_BINS" -eq 0 ]; then
    echo "ERROR: No *.fa files found in $BINS_DIR/"
    exit 1
fi

NUM_CHUNK=$(( NUM_BINS / 3 ))
(( NUM_CHUNK < 1 )) && NUM_CHUNK=1      # If NUM_BINS < 3
echo "Detected $NUM_BINS bins → chunk size = $NUM_CHUNK"


DREP_GENOMES_DIR="$DREP_OUT/dereplicated_genomes"

if [[ -d "$DREP_GENOMES_DIR" && $(ls -1 "$DREP_GENOMES_DIR"/*.fa 2>/dev/null | wc -l) -gt 0 ]]; then
    echo "dRep output already exists → skipping dereplication."
else
    echo "Running dRep dereplicate..."
    dRep check_dependencies || echo "WARNING: Missing optional dependencies, continuing..." # not using -g "$BINS_DIR"/*.fa
    dRep dereplicate "$DREP_OUT" \
        -g "$GENOME_LIST" \
        -p 56 \
        -l 10000 \
        -pa 0.9 \
        -sa 0.95 \
        -nc 0.5 \
        -cm larger \
        --S_algorithm fastANI \
        -comp 50 \
        --primary_chunksize "$NUM_CHUNK" \
        --run_tertiary_clustering
fi

# Verify dRep output exists and is non-empty
#NUM_MAGs=$(ls -1 "$DREP_OUT/dereplicated_genomes"/*.fa 2>/dev/null | wc -l) 
MAG_LIST="$DREP_OUT/mag_list.txt"
find -L "$DREP_OUT/dereplicated_genomes" -maxdepth 1 -type f -name "*.fa" > "$MAG_LIST"
NUM_MAGs=$(wc -l < "$MAG_LIST")

if [[ "$NUM_MAGs" -eq 0 ]]; then
    echo "ERROR: No dereplicated MAGs produced."
    exit 1
fi

echo "dRep produced $NUM_MAGs MAGs."


echo "Running gtdbtk_r226..."
#export GTDBTK_DATA_PATH=/data/MicrobiomeCore/databases/gtdbtk_r226/gtdbtk_r226_data/release226
#echo $GTDBTK_DATA_PATH
# add condition where the output folder is changed to gtdbtk_custom_out if default ref db is not used and provided by data path

if [[ -d "$DREP_OUT/gtdbtk_r226_out" ]]; then
    echo "GTDB-Tk output exists → skipping classification."
else
    mkdir -p "$DREP_OUT/gtdbtk_r226_out"

    gtdbtk classify_wf \
  --genome_dir "$DREP_OUT/dereplicated_genomes" \
  --out_dir "$DREP_OUT/gtdbtk_r226_out" \
  --skip_ani_screen \
  --cpus 32 \
  --pplacer_cpus 8 -x fa &> "$DREP_OUT/gtdbtk_r226_out/gtdbtk.log"
fi


echo "DONE."
