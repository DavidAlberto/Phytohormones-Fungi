#!/usr/bin/env bash
# =============================================================================
# Script:      extract_hmm_ids.sh
# Description: Fetches individual HMM profiles from a local Pfam-A database
#              for a list of Pfam accessions. One .hmm file is generated per
#              accession in the output directory.
#
# Usage:       bash extract_hmm_ids.sh <Pfam-A.hmm> <pfam_list.txt> <output_dir>
#
#   <Pfam-A.hmm>     — path to local Pfam-A HMM database
#   <pfam_list.txt>  — text file with one Pfam accession per line (e.g. PF00117)
#   <output_dir>     — directory where individual .hmm files will be written
#
# Output:      <output_dir>/<PFAM_ID>.hmm   — one HMM profile per accession
#              <output_dir>/error_log.txt   — accessions not found in database
#
# Dependencies: hmmfetch (HMMER 3.x)
# =============================================================================

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <Pfam-A.hmm> <pfam_list.txt> <output_dir>"
    exit 1
fi

HMM_DB="$1"
LIST="$2"
OUTPUT_DIR="$3"
LOG_FILE="${OUTPUT_DIR}/error_log.txt"

if ! command -v hmmfetch &>/dev/null; then
    echo "Error: hmmfetch not found. Ensure HMMER is installed and in PATH."
    exit 1
fi

if [ ! -f "$HMM_DB" ]; then
    echo "Error: Pfam HMM database not found: '$HMM_DB'"
    exit 1
fi

if [ ! -f "$LIST" ]; then
    echo "Error: Pfam accession list not found: '$LIST'"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

SUCCESS_COUNT=0
ERROR_COUNT=0

while read -r PFAM; do
    [[ -z "$PFAM" || "$PFAM" =~ ^# ]] && continue
    echo "Fetching $PFAM..."
    if hmmfetch "$HMM_DB" "$PFAM" > "${OUTPUT_DIR}/${PFAM}.hmm" 2>/dev/null; then
        ((SUCCESS_COUNT++))
    else
        echo "Error: $PFAM not found in $HMM_DB" | tee -a "$LOG_FILE"
        ((ERROR_COUNT++))
    fi
done < "$LIST"

echo "Done."
echo "  ✅ Fetched successfully : $SUCCESS_COUNT"
echo "  ❌ Not found (errors)   : $ERROR_COUNT"
