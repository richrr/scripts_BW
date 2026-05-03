#!/usr/bin/env bash
set -euo pipefail

# /bin/ls *.fastq.gz | sed 's/_1.fastq.gz//; s/_2.fastq.gz//' | sort -u > ../completed_samples.txt
# bash /data/rodriguesrr/scripts/bash/build_download_swarm_missed_fetchngs.sh metadata fastq completed_samples.txt
#Run swarm:
#  swarm -f download_fastq.swarm -g 8 -t 4 --time 4:00:00 --logdir logs
#Check:
#  md5sum -c checksums.txt


METADATA_DIR=$1
OUTDIR=$2
COMPLETED_FILE=${3:-""}
SWARM_FILE=${4:-download_fastq.swarm}
MD5_FILE=${5:-checksums.txt}

mkdir -p "$OUTDIR"

> "$SWARM_FILE"
> "$MD5_FILE"

echo "Metadata dir:   $METADATA_DIR"
echo "Output dir:     $OUTDIR"
echo "Completed file: ${COMPLETED_FILE:-NONE}"
echo "Swarm file:     $SWARM_FILE"
echo "MD5 file:       $MD5_FILE"

for f in "$METADATA_DIR"/*.tsv; do
  awk -F'\t' -v outdir="$OUTDIR" \
              -v swarm="$SWARM_FILE" \
              -v md5file="$MD5_FILE" \
              -v donefile="$COMPLETED_FILE" '
  BEGIN {OFS="\t"}

  NR==1 {
    for (i=1;i<=NF;i++) {
      if ($i=="id") idcol=i
      if ($i=="fastq_1") fq1col=i
      if ($i=="fastq_2") fq2col=i
      if ($i=="md5_1") md5_1col=i
      if ($i=="md5_2") md5_2col=i
    }
    next
  }

  {
    id = $(idcol)
    fq1 = $(fq1col)
    fq2 = $(fq2col)
    m1  = $(md5_1col)
    m2  = $(md5_2col)

    if (fq1=="" || fq2=="") next

    file1 = id "_1.fastq.gz"
    file2 = id "_2.fastq.gz"

    # check if already completed (if file provided)
    completed = 0
    if (donefile != "") {
      cmd = "grep -qx \"" id "\" " donefile
      if (system(cmd) == 0) completed = 1
    }

    # ALWAYS write checksum entries
    print m1 "  " outdir "/" file1 >> md5file
    print m2 "  " outdir "/" file2 >> md5file

    # only queue download if not completed
    if (!completed) {
      print "wget -c -t 3 -T 120 " fq1 " -O " outdir "/" file1 >> swarm
      print "wget -c -t 3 -T 120 " fq2 " -O " outdir "/" file2 >> swarm
    }
  }
  ' "$f"
done

echo "Done."
echo "Run swarm:"
echo "  swarm -f $SWARM_FILE -g 8 -t 4 --time 4:00:00 --logdir logs"
echo "Check:"
echo "  md5sum -c $MD5_FILE"