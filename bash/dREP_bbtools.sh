# this script has some hardcoded parts assuming the reads were cleaned using mag pipiline and contain "_phix_removed.unmapped_?.fastq.gz"
# it also does not have gtdb-tk on the dRep MAGs.
# -- use the script inside nfcore.mag.accessories --

# ml dRep
# dRep check_dependencies
# NUM_BINS=$(ls 2>/dev/null -Ubad1 -- input.directory.bins/*fa | wc -l)
# NUM_CHUNK=$(( $NUM_BINS/3 ))
# dRep dereplicate -d -g <<input.directory.bins/*fa>> -p 56 -l 10000 -pa 0.9 -sa 0.95 -nc 0.5 -cm 'larger' --S_algorithm fastANI -comp 50 --primary_chunksize $NUM_CHUNK --run_tertiary_clustering <<output.directory.path>>



# module load bbtools/39.06

# # indexing
# unset _JAVA_OPTIONS
# bbtools bbsplit -Xmx64g ref=GenomeBinning/DASTool/bins.dRep.MAGs/dereplicated_genomes/ path=GenomeBinning/DASTool/bins.dRep.MAGs/dereplicated_genomes.ref


# #quantification

# bbtools bbsplit -Xmx64g unpigz=t threads=8 minid=0.90 sortscafs=f nzo=f \
# ref=GenomeBinning/DASTool/bins.dRep.MAGs/dereplicated_genomes/ path=GenomeBinning/DASTool/bins.dRep.MAGs/dereplicated_genomes.ref \
# statsfile=GenomeBinning/DASTool/bins.dRep.samples/MF0015pGTSH10d20221212.statsfile \
# scafstats=GenomeBinning/DASTool/bins.dRep.samples/MF0015pGTSH10d20221212.scafstats \
# covstats=GenomeBinning/DASTool/bins.dRep.samples/MF0015pGTSH10d20221212.covstat \
# rpkm=GenomeBinning/DASTool/bins.dRep.samples/MF0015pGTSH10d20221212.rpkm \
# refstats=GenomeBinning/DASTool/bins.dRep.samples/MF0015pGTSH10d20221212.refstats \
# in=QC_shortreads/remove_phix/MF0015pGTSH10d20221212_run0_phix_removed.unmapped_1.fastq.gz \
# in2=QC_shortreads/remove_phix/MF0015pGTSH10d20221212_run0_phix_removed.unmapped_2.fastq.gz

# cd /data/rodriguesrr/MetaflowX/test_pitts_mag_no_spades_busco_metabinner

# bash /data/rodriguesrr/scripts/bash/dREP_bbtools.sh GenomeBinning/DASTool/bins ../test.pitts.samplesheet.mag.csv QC_shortreads/remove_phix/

# cd /data/rodriguesrr/Koltsova/analysis/Nov2025_IL22_Alb_Vil/nf-core-mag/odir_mag_no_spades_busco_metabinner
#swarm --module=drep/3.4.2,bbtools/39.06 -g 300 -t 24 --time=2-00:00:00 --logdir swarmlogs -f run_dREP_bbtools.swarm


#!/usr/bin/env bash

set -euo pipefail

# ---------------------------
# USAGE
# ---------------------------
if [ $# -ne 3 ]; then
    echo "Usage: $0 <bins_directory> <sample_sheet> <clean_fastq_directory>"
    exit 1
fi

BINS_DIR=$1            # e.g., GenomeBinning/DASTool/bins
SAMPLE_SHEET=$2        # file with sample names, first column used
FASTQ_DIR=$3           # e.g., QC_shortreads/remove_phix

# ---------------------------
# OUTPUT DIRECTORIES
# ---------------------------
DREP_OUT="${BINS_DIR}.dRep.MAGs"
SAMPLE_OUT="${BINS_DIR}.dRep.samples"

mkdir -p "$DREP_OUT" "$SAMPLE_OUT"

echo "Bins directory:           $BINS_DIR"
echo "Dereplicated bins out:    $DREP_OUT"
echo "Sample output directory:  $SAMPLE_OUT"
echo "FASTQ directory:          $FASTQ_DIR"
echo "Sample sheet:             $SAMPLE_SHEET"

# ---------------------------
# LOAD MODULES
# ---------------------------
#module load sqlite/3.51.1

# tell Lmod not to fail the script if a module complains
#set +e
module load drep/3.4.2
#set -e

module load bbtools/39.06

# ---------------------------
# DEREPLICATION
# ---------------------------
echo "Counting bins..."
NUM_BINS=$(ls -1 "$BINS_DIR"/*.fa 2>/dev/null | wc -l)

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
    dRep check_dependencies || echo "WARNING: Missing optional dependencies, continuing..."
    dRep dereplicate "$DREP_OUT" \
        -g "$BINS_DIR"/*.fa \
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
NUM_MAGs=$(ls -1 "$DREP_OUT/dereplicated_genomes"/*.fa 2>/dev/null | wc -l)

if [[ "$NUM_MAGs" -eq 0 ]]; then
    echo "ERROR: No dereplicated MAGs produced."
    exit 1
fi

echo "dRep produced $NUM_MAGs MAGs."


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
# PER-SAMPLE QUANTIFICATION
# ---------------------------
echo "Starting per-sample quantification..."

# Read CSV, skip header, extract first two columns
while IFS=',' read -r SAMPLE RUN _; do
    # Match files like:
    # MF0015pGTSH10d20221212_run0_phix_removed.unmapped_1.fastq.gz
    R1="${FASTQ_DIR}/${SAMPLE}_run${RUN}_phix_removed.unmapped_1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_run${RUN}_phix_removed.unmapped_2.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "WARNING: FASTQs not found for $SAMPLE (run $RUN) → skipping."
        continue
    fi

    RPKM_OUT="$SAMPLE_OUT/${SAMPLE}.rpkm"
    if [[ -s "$RPKM_OUT" ]]; then
        echo "Sample $SAMPLE already quantified → skipping."
        continue
    fi

    
    echo "Processing sample: $SAMPLE"
    echo "  R1: $R1"
    echo "  R2: $R2"

    unset _JAVA_OPTIONS
    LOG="$SAMPLE_OUT/${SAMPLE}.bbsplit.log"

    bbtools bbsplit -Xmx64g \
      unpigz=t \
      threads=8 \
      minid=0.90 \
      sortscafs=f \
      nzo=f \
      ref="$DREP_OUT/dereplicated_genomes/" \
      path="$DREP_OUT/dereplicated_genomes.ref" \
      statsfile="$SAMPLE_OUT/${SAMPLE}.statsfile" \
      scafstats="$SAMPLE_OUT/${SAMPLE}.scafstats" \
      covstats="$SAMPLE_OUT/${SAMPLE}.covstat" \
      rpkm="$SAMPLE_OUT/${SAMPLE}.rpkm" \
      refstats="$SAMPLE_OUT/${SAMPLE}.refstats" \
      in="$R1" \
      in2="$R2" &> "$LOG"

done < <(tail -n +2 "$SAMPLE_SHEET"; echo)

# Verify per-sample outputs at the very end
EXPECTED=$(tail -n +2 "$SAMPLE_SHEET" | wc -l)
OBSERVED=$(ls -1 "$SAMPLE_OUT"/*.rpkm 2>/dev/null | wc -l)
echo "Expected samples: $EXPECTED"
echo "Observed RPKM files: $OBSERVED"

if [[ "$OBSERVED" -ne "$EXPECTED" ]]; then
    echo "WARNING: Not all samples were successfully quantified."
    echo "You can safely rerun the job to resume."
else
    echo "All samples successfully quantified."
fi


echo "DONE."
