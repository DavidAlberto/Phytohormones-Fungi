#!/usr/bin/env bash
# =============================================================================
# Script:      search_hmm_ids.sh
# Description: Searches all HMM profiles in a directory against a target
#              proteome using hmmsearch. Produces one domain hit table per
#              profile in the --domtblout format for downstream domain
#              coordinate extraction.
#
# Usage:       bash search_hmm_ids.sh <hmm_dir> <proteome.fasta> <output_dir>
#
#   <hmm_dir>        — directory containing .hmm files
#   <proteome.fasta> — protein FASTA to search (e.g., T. atroviride v3)
#   <output_dir>     — directory where hit tables will be written
#
# Output:      <output_dir>/<PFAM_ID>_hmmer.tbl — domain hit table per profile
#
# Dependencies: hmmsearch (HMMER 3.x)
# =============================================================================

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <hmm_dir> <proteome.fasta> <output_dir>"
    exit 1
fi

HMM_FOLDER="$1"
PROTEOME="$2"
OUTPUT_DIR="$3"

if ! command -v hmmsearch &>/dev/null; then
    echo "Error: hmmsearch not found. Ensure HMMER is installed and in PATH."
    exit 1
fi

if [ ! -d "$HMM_FOLDER" ]; then
    echo "Error: HMM directory not found: '$HMM_FOLDER'"
    exit 1
fi

if [ ! -f "$PROTEOME" ]; then
    echo "Error: Proteome file not found: '$PROTEOME'"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

for HMM_FILE in "$HMM_FOLDER"/*.hmm; do
    BASENAME=$(basename "$HMM_FILE" .hmm)
    OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}_hmmer.tbl"
    echo "Searching with $BASENAME..."
    hmmsearch --domtblout "$OUTPUT_FILE" "$HMM_FILE" "$PROTEOME"
done

echo "Search complete. Results saved in '${OUTPUT_DIR}'."
