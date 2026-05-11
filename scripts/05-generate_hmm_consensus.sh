#!/usr/bin/env bash
# =============================================================================
# Script:      generate_hmm_consensus.sh
# Description: Generates consensus sequences from HMM profiles using hmmemit.
#              Processes all .hmm files in the input directory and writes one
#              consensus FASTA per profile to the output directory.
#
# Usage:       bash generate_hmm_consensus.sh <hmm_dir> <output_dir>
#
#   <hmm_dir>    — directory containing .hmm files (output of extract_hmm_ids.sh)
#   <output_dir> — directory where consensus FASTA files will be written
#
# Output:      <output_dir>/<PFAM_ID>_consensus.fasta — one file per profile
#
# Dependencies: hmmemit (HMMER 3.x)
# =============================================================================

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <hmm_dir> <output_dir>"
    exit 1
fi

HMM_FOLDER="$1"
OUTPUT_DIR="$2"

if ! command -v hmmemit &>/dev/null; then
    echo "Error: hmmemit not found. Ensure HMMER is installed and in PATH."
    exit 1
fi

if [ ! -d "$HMM_FOLDER" ]; then
    echo "Error: HMM directory not found: '$HMM_FOLDER'"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

for HMM_FILE in "$HMM_FOLDER"/*.hmm; do
    BASENAME=$(basename "$HMM_FILE" .hmm)
    OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}_consensus.fasta"
    echo "Generating consensus for $BASENAME..."
    hmmemit -c "$HMM_FILE" > "$OUTPUT_FILE"
done

echo "Done. Consensus sequences saved in '${OUTPUT_DIR}'."
