#!/usr/bin/env bash
# =============================================================================
# Script:      prepare_blastdb.sh
# Description: Standardizes FASTA headers from downloaded proteomes and merges
#              them into a single non-redundant database ready for BLASTp.
#
#              Each sequence header is prefixed with the organism's TaxID
#              (extracted from the filename), producing the format:
#                  >TaxID|original_header
#              This prefix is required for organism extraction in downstream
#              filtering (Blastp-analysis.Rmd, ORGANISM_METHOD = "taxid").
#
# Usage:       bash prepare_blastdb.sh
#              (edit INPUT_DIR / OUTPUT_DIR / FINAL_DB_FILE below if needed)
#
# Input:       proteomes/<TaxID>.fasta  — one file per organism, named by TaxID
# Output:      proteomes_modified_headers/<TaxID>.fasta  — renamed headers
#              hormoneDB.fasta  — merged multifasta database
#
# Next step:   makeblastdb -in hormoneDB.fasta -dbtype prot -parse_seqids \
#                          -out blastDB
#
# Dependencies: bash, awk
# =============================================================================

set -euo pipefail

# --- Configuration ---
readonly DATA_DIR="data"
readonly INPUT_DIR="${DATA_DIR}/proteomes"
readonly OUTPUT_DIR="${DATA_DIR}/proteomes_modified_headers"
readonly FINAL_DB_FILE="${DATA_DIR}/hormoneDB.fasta"

# --- Colored output (suppressed when stdout is not a TTY) ---
if [ -t 1 ]; then
    readonly GREEN='\033[0;32m'
    readonly BLUE='\033[0;34m'
    readonly NC='\033[0m'
else
    readonly GREEN='' BLUE='' NC=''
fi

log_info()    { echo -e "${BLUE}[INFO]${NC} $*"; }
log_success() { echo -e "${GREEN}[OK]${NC} $*"; }

# --- Main ---
log_info "Starting BLAST database preparation."
log_info "Input directory : ${INPUT_DIR}"
log_info "Output directory: ${OUTPUT_DIR}"

if [ ! -d "${INPUT_DIR}" ]; then
    echo "[ERROR] Input directory '${INPUT_DIR}' not found." >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# Process each FASTA file: prefix every sequence header with the TaxID
# (derived from the filename, e.g., 5476.fasta → TaxID 5476).
for fasta_file in "${INPUT_DIR}"/*.fasta; do
    [ -f "$fasta_file" ] || continue

    FILENAME=$(basename "$fasta_file")
    TAXID="${FILENAME%.fasta}"
    OUTPUT_FILE="${OUTPUT_DIR}/${FILENAME}"

    echo "  Processing: ${FILENAME}  →  adding TaxID prefix: ${TAXID}"

    awk -v taxid="$TAXID" \
        '/^>/ { sub(/^>/, ">" taxid "|"); print; next } { print }' \
        "$fasta_file" > "$OUTPUT_FILE"
done

log_success "All files processed and saved to '${OUTPUT_DIR}'."

# Merge all renamed proteomes into a single database file
log_info "Merging all proteomes into: ${FINAL_DB_FILE}"
cat "${OUTPUT_DIR}"/*.fasta > "${FINAL_DB_FILE}"

log_success "Database file ready: ${FINAL_DB_FILE}"
log_info  "Create the BLAST index with:"
echo "    makeblastdb -in ${FINAL_DB_FILE} -dbtype prot -parse_seqids"
