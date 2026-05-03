#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <metadata_dir> <output_dir> [num_threads]"
  exit 1
fi

METADATA_DIR="$1"
OUTDIR="$2"
THREADS="${3:-4}"

ASPERA_KEY="/opt/aspera/aspera_tokenauth_id_rsa"

mkdir -p "$OUTDIR"

echo "Metadata dir: $METADATA_DIR"
echo "Output dir:   $OUTDIR"
echo "Threads:      $THREADS"
echo "Using key:    $ASPERA_KEY"
echo

TMP_LINKS=$(mktemp)
TMP_MD5=$(mktemp)

# ============================
# Parse metadata
# ============================
for f in "$METADATA_DIR"/*.tsv; do
  awk -F'\t' '
    NR==1 {
      for (i=1; i<=NF; i++) {
        if ($i=="fastq_aspera") col_aspera=i
        if ($i=="fastq_ftp") col_ftp=i
        if ($i=="md5_1") col_md5_1=i
        if ($i=="md5_2") col_md5_2=i
      }
    }
    NR==2 {
      # Aspera links
      if (col_aspera && $col_aspera != "" && $col_aspera != "NA") {
        print $col_aspera >> "'$TMP_LINKS'"
      }

      # MD5
      n = split($col_ftp, ftp_arr, ";")

      if (n >= 1 && col_md5_1 && $col_md5_1 != "") {
        split(ftp_arr[1], a, "/")
        print $col_md5_1 "  " a[length(a)] >> "'$TMP_MD5'"
      }

      if (n == 2 && col_md5_2 && $col_md5_2 != "") {
        split(ftp_arr[2], a, "/")
        print $col_md5_2 "  " a[length(a)] >> "'$TMP_MD5'"
      }
    }
  ' "$f"
done

ASPERA_LIST=$(mktemp)
cat "$TMP_LINKS" | tr ';' '\n' | sed 's#^#era-fasp@#' > "$ASPERA_LIST"

TOTAL=$(wc -l < "$ASPERA_LIST")
echo "Total FASTQ files to process: $TOTAL"

CHECKSUM_FILE="$OUTDIR/checksums.txt"
sort -u "$TMP_MD5" > "$CHECKSUM_FILE"

echo "Checksum file written to: $CHECKSUM_FILE"
echo

# ============================
# Run Aspera
# ============================
cd "$OUTDIR"

cat "$ASPERA_LIST" | xargs -I {} -P "$THREADS" bash -c '
  echo "Downloading: {}"
  ascp -QT -l 500M -P33001 -k2 -i "'"$ASPERA_KEY"'" {} .
'

echo
echo "Download complete."

#rm -f "$TMP_LINKS" "$TMP_MD5" "$ASPERA_LIST"

echo
echo "To verify downloads:"
echo "cd $OUTDIR && md5sum -c checksums.txt"