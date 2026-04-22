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

# ---------------------------
# DEREPLICATION INPUT RESOLUTION
# ---------------------------

PARENT_DIR="$(dirname "$BINS_DIR")"

GENOME_LIST="$PARENT_DIR/genomes_list.txt"
CHECKM2_FILTERED="$PARENT_DIR/checkm2_filtered.tsv"


echo "Resolving genome inputs for dRep..."

# ---------------------------
# CASE 1: Already aggregated (preferred path)
# ---------------------------
if [[ -f "$GENOME_LIST" && -s "$GENOME_LIST" && \
      -f "$CHECKM2_FILTERED" && -s "$CHECKM2_FILTERED" ]]; then
    echo "Detected aggregated run (genomes_list + filtered QC exist)."

# ---------------------------
# CASE 2: Single batch or missing gather outputs
# ---------------------------
else
    echo "Aggregated files missing → building QC-filtered genome list from CheckM2..."

    #CHECKM2_SUMMARY="GenomeBinning/QC/checkm2_summary.tsv"
    CHECKM2_SUMMARY="$(dirname "$BINS_DIR")/../QC/checkm2_summary.tsv"
    #[[ -f "$CHECKM2_SUMMARY" ]] || { echo "ERROR: CheckM2 not found"; exit 1; }


    # sanity check input exists
    if [[ ! -s "$CHECKM2_SUMMARY" ]]; then
        echo "ERROR: CheckM2 summary not found at $CHECKM2_SUMMARY"
        exit 1
    fi


    HEADER="$CHECKM2_SUMMARY"

    COMP_COL=$(head -1 "$HEADER" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="Completeness") print i}')
    CONT_COL=$(head -1 "$HEADER" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="Contamination") print i}')
    NAME_COL=$(head -1 "$HEADER" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="Name") print i}')

    if [[ -z "$COMP_COL" || -z "$CONT_COL" || -z "$NAME_COL" ]]; then
        echo "ERROR: Missing required columns in CheckM2 file"
        exit 1
    fi
    
    # ---------------------------
    # build QC filtered table
    # (keep your thresholds consistent with pipeline policy)
    # ---------------------------
    awk -F'\t' -v comp="$COMP_COL" -v cont="$CONT_COL" '
NR==1 || ($comp >= 50 && $cont <= 10)
' "$CHECKM2_SUMMARY" > "$CHECKM2_FILTERED"

    # ---------------------------
    # build genome list from QC table
    # ---------------------------
    awk -F'\t' -v name="$NAME_COL" -v dir="$BINS_DIR" '
NR>1 {print dir "/" $name ".fa"}
' "$CHECKM2_FILTERED" > "$GENOME_LIST"

fi

# ---------------------------
# FINAL VALIDATION
# ---------------------------
NUM_BINS=$(wc -l < "$GENOME_LIST")

if [[ "$NUM_BINS" -eq 0 ]]; then
    echo "ERROR: Genome list is empty after resolution step."
    exit 1
fi

echo "Final genome count for dRep: $NUM_BINS"



if [[ $(wc -l < "$GENOME_LIST") -ne $(awk 'NR>1' "$CHECKM2_FILTERED" | wc -l) ]]; then
    echo "WARNING: mismatch between QC table and genome list"
fi




# dynamic chunk sizing (keep your logic but safer scaling)
#NUM_CHUNK=$(( NUM_BINS / 3 ))
#(( NUM_CHUNK < 1 )) && NUM_CHUNK=1
NUM_CHUNK=$(( NUM_BINS / 10 ))
(( NUM_CHUNK < 1000 )) && NUM_CHUNK=1000
(( NUM_CHUNK > 5000 )) && NUM_CHUNK=5000
echo "Chunk size = $NUM_CHUNK"



#DREP_GENOMES_DIR="$DREP_OUT/dereplicated_genomes"
#if [[ -d "$DREP_GENOMES_DIR" && $(ls -1 "$DREP_GENOMES_DIR"/*.fa 2>/dev/null | wc -l) -gt 0 ]]; then
if [[ -s "$DREP_OUT/data_tables/Cdb.csv" && \
      -d "$DREP_OUT/dereplicated_genomes" && \
      $(find "$DREP_OUT/dereplicated_genomes" -name "*.fa" | wc -l) -gt 0 ]]; then
    echo "dRep output already exists → skipping dereplication."
else
    echo "Running dRep dereplicate..."
    dRep check_dependencies || echo "WARNING: Missing optional dependencies, continuing..." # not using -g "$BINS_DIR"/*.fa

    LARGE_THRESHOLD=15000   # switch point (tunable)

    if [[ "$NUM_BINS" -lt "$LARGE_THRESHOLD" ]]; then
        echo "Small/medium dataset → running single-pass dRep"

        # replacing -comp 50 with --ignoreGenomeQuality
        dRep dereplicate "$DREP_OUT" \
            -g "$GENOME_LIST" \
            -p 56 \
            -l 10000 \
            -pa 0.9 \
            -sa 0.95 \
            -nc 0.5 \
            -cm larger \
            --S_algorithm fastANI \
            --ignoreGenomeQuality \
            --primary_chunksize "$NUM_CHUNK" \
            --run_tertiary_clustering
    else
        echo "Large dataset detected ($NUM_BINS genomes) → using 2-stage dRep"

        STAGE1_DIR="$DREP_OUT/stage1_chunks"
        mkdir -p "$STAGE1_DIR"

        echo "Splitting genome list..."
        split -l 5000 "$GENOME_LIST" "$STAGE1_DIR/chunk_"

        echo "Running Stage 1 dRep on chunks..."

        for chunk in "$STAGE1_DIR"/chunk_*; do
            chunk_name=$(basename "$chunk")
            outdir="$STAGE1_DIR/${chunk_name}.drep"

            echo "Processing $chunk_name"

            dRep dereplicate "$outdir" \
                -g "$chunk" \
                -p 16 \
                -pa 0.9 \
                -sa 0.9 \
                --S_algorithm fastANI \
                --ignoreGenomeQuality
        done

        echo "Collecting Stage 1 representatives..."

        STAGE2_INPUT="$DREP_OUT/stage2_genomes.txt"
        find "$STAGE1_DIR" -path "*/dereplicated_genomes/*.fa" > "$STAGE2_INPUT"

        echo "Running Stage 2 dRep..."

        dRep dereplicate "$DREP_OUT" \
            -g "$STAGE2_INPUT" \
            -p 32 \
            -l 10000 \
            -pa 0.9 \
            -sa 0.95 \
            -nc 0.5 \
            -cm larger \
            --S_algorithm fastANI \
            --ignoreGenomeQuality \
            --primary_chunksize 1000 \
            --run_tertiary_clustering
    fi
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
