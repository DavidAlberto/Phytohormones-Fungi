#!/usr/bin/env bash
# =============================================================================
# Script:      download_proteomes.sh
# Description: Downloads proteomes for all organisms listed in organisms.tsv
#              using a hierarchical fallback strategy to maximize annotation
#              quality. Each organism is attempted in this order:
#                1. UniProt reference proteome
#                2. UniProt reviewed sequences (Swiss-Prot)
#                3. UniProt complete sequences (TrEMBL)
#                4. NCBI RefSeq (fallback if all UniProt attempts fail)
#
# Usage:       bash download_proteomes.sh
#
# Input:       organisms.tsv                  — one scientific name per line
# Output:      proteomes/<TaxID>.fasta        — one file per organism
#              manifests/taxids.tsv           — TaxID resolution log
#              manifests/proteomes.tsv        — download manifest
#              manifests/errors.log           — errors and warnings
#              manifests/statistics.txt       — final summary
#
# Dependencies: curl, esearch, efetch (NCBI E-Direct), awk, jq (optional)
#
# Notes:
#   - Idempotent: skips organisms whose FASTA already exists in proteomes/.
#   - T. atroviride v3 (unpublished) must be added manually after this script.
#   - jq improves JSON parsing reliability; the script falls back to grep
#     if jq is not available.
# =============================================================================

# Do NOT use set -e: we want to continue even if individual organisms fail.
set -uo pipefail

# --- Configuration ---
readonly DATA_DIR="data"
readonly ORGANISM_LIST="${DATA_DIR}/organisms.tsv"
readonly MANIFEST_DIR="${DATA_DIR}manifests"
readonly PROTEOME_DIR="${DATA_DIR}proteomes"
readonly TAXID_MANIFEST="${MANIFEST_DIR}/taxids.tsv"
readonly PROTEOME_MANIFEST="${MANIFEST_DIR}/proteomes.tsv"
readonly ERROR_LOG="${MANIFEST_DIR}/errors.log"
readonly STATS_FILE="${MANIFEST_DIR}/statistics.txt"

readonly MAX_RETRIES=3
readonly TIMEOUT=300

# --- Colored output ---
if [ -t 1 ]; then
    readonly RED='\033[0;31m'
    readonly GREEN='\033[0;32m'
    readonly YELLOW='\033[1;33m'
    readonly BLUE='\033[0;34m'
    readonly NC='\033[0m'
else
    readonly RED='' GREEN='' YELLOW='' BLUE='' NC=''
fi

log_info()    { echo -e "${BLUE}[INFO]${NC} $*"; }
log_success() { echo -e "${GREEN}[OK]${NC} $*"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $*" >&2
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >> "$ERROR_LOG"; }
log_warn()    { echo -e "${YELLOW}[WARN]${NC} $*"
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN: $*" >> "$ERROR_LOG"; }

# --- Initialization ---
init_environment() {
    mkdir -p "${MANIFEST_DIR}" "${PROTEOME_DIR}"
    echo -e "Organism\tStatus\tDetails" > "${TAXID_MANIFEST}"
    echo -e "Organism\tTaxID\tProteomeID\tType\tSeqCount\tFilename\tStatus\tSource"  > "${PROTEOME_MANIFEST}"
    > "$ERROR_LOG"

    echo "========================================================================"
    log_info "Bulk proteome download — UniProt with NCBI fallback"
    log_info "Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "========================================================================"
    log_info "Output locations:"
    log_info "  Proteomes      : ${PROTEOME_DIR}/"
    log_info "  TaxID manifest : ${TAXID_MANIFEST}"
    log_info "  Proteome log   : ${PROTEOME_MANIFEST}"
    log_info "  Error log      : ${ERROR_LOG}"
    echo "------------------------------------------------------------------------"
}

check_dependencies() {
    local missing=0
    for cmd in esearch efetch curl awk; do
        if ! command -v "$cmd" &> /dev/null; then
            log_error "Missing dependency: $cmd"
            missing=1
        fi
    done
    if ! command -v jq &> /dev/null; then
        log_warn "jq not found (optional — improves JSON parsing reliability)"
    fi
    if [ $missing -eq 1 ]; then
        log_error "Install missing dependencies:"
        echo "  Ubuntu/Debian: sudo apt install ncbi-entrez-direct curl gawk jq"
        exit 1
    fi
}

# --- TaxID resolution ---
# Queries NCBI Taxonomy by scientific name; retries up to MAX_RETRIES times.
get_taxid() {
    local organism="$1"
    local attempt=1
    local taxid_results=""
    local taxid_count

    while [ $attempt -le $MAX_RETRIES ]; do
        taxid_results=$(esearch -db taxonomy -query "${organism}[SCIN]" 2>/dev/null | \
                        efetch -format uid 2>/dev/null | \
                        grep -E '^[0-9]+$' 2>/dev/null || true)
        [ -n "$taxid_results" ] && break
        log_warn "Attempt $attempt/$MAX_RETRIES failed for: $organism"
        attempt=$((attempt + 1))
        sleep 2
    done

    taxid_count=$(echo "$taxid_results" | grep -c -E "^[0-9]*$")

    if [ "$taxid_count" -eq 1 ]; then
        echo "$taxid_results"
        echo -e "${organism}\t${taxid_results}\tOK" >> "${TAXID_MANIFEST}"
        return 0
    elif [ "$taxid_count" -eq 0 ]; then
        log_error "TaxID not found: $organism"
        echo -e "${organism}\tNOT_FOUND\t-" >> "${TAXID_MANIFEST}"
        return 1
    else
        log_warn "Ambiguous name: $organism (${taxid_count} matches)"
        echo -e "${organism}\tAMBIGUOUS\t${taxid_count}_matches" >> "${TAXID_MANIFEST}"
        return 1
    fi
}

# --- UniProt download ---
# Attempts reference proteome → reviewed → complete for the given TaxID.
get_reference_proteome() {
    local taxid="$1"
    local temp_json="${PROTEOME_DIR}/.temp_${taxid}.json"

    curl -s -G "https://rest.uniprot.org/proteomes/search" \
        --data-urlencode "query=organism_id:${taxid} AND reference:true" \
        --data-urlencode "format=json" \
        --max-time 30 \
        -o "$temp_json" 2>/dev/null

    local proteome_id=""
    if command -v jq &> /dev/null; then
        proteome_id=$(jq -r '.results[0].id // empty' "$temp_json" 2>/dev/null || echo "")
    else
        proteome_id=$(grep -oP '"id"\s*:\s*"\K[^"]+' "$temp_json" 2>/dev/null | head -1 || echo "")
    fi

    rm -f "$temp_json"
    echo "$proteome_id"
}

download_proteome_uniprot() {
    local organism="$1"
    local taxid="$2"
    local out_fa="${PROTEOME_DIR}/${taxid}.fasta"

    # Skip if already downloaded
    if [ -f "${out_fa}" ]; then
        log_info "  ✓ Already exists: ${out_fa}. Skipping."
        echo -e "${organism}\t${taxid}\t-\tSKIPPED_EXISTS\t-\t${out_fa}\tSUCCESS" >> "${PROTEOME_MANIFEST}"
        return 0
    fi

    local out_tmp="${out_fa}.tmp"
    local proteome_type="UNKNOWN"
    local proteome_id="-"
    local seq_count=0

    log_info "  Searching proteome for TaxID: ${taxid}"

    # Strategy 1: reference proteome
    proteome_id=$(get_reference_proteome "$taxid")
    if [ -n "$proteome_id" ]; then
        log_info "  → Downloading REFERENCE proteome: ${proteome_id}"
        curl -s -G "https://rest.uniprot.org/uniprotkb/stream" \
            --data-urlencode "query=proteome:${proteome_id}" \
            --data-urlencode "format=fasta" \
            --compressed --max-time $TIMEOUT \
            -o "${out_tmp}" 2>>"$ERROR_LOG"
        grep -q "^>" "${out_tmp}" 2>/dev/null && proteome_type="REFERENCE"
    fi

    # Strategy 2: reviewed sequences (Swiss-Prot)
    if [ ! -s "${out_tmp}" ] || ! grep -q "^>" "${out_tmp}"; then
        log_info "  → Downloading REVIEWED sequences (Swiss-Prot)"
        curl -s -G "https://rest.uniprot.org/uniprotkb/stream" \
            --data-urlencode "query=organism_id:${taxid} AND reviewed:true" \
            --data-urlencode "format=fasta" \
            --compressed --max-time $TIMEOUT \
            -o "${out_tmp}" 2>>"$ERROR_LOG"
        grep -q "^>" "${out_tmp}" 2>/dev/null && proteome_type="REVIEWED"
    fi

    # Strategy 3: complete proteome (TrEMBL)
    if [ ! -s "${out_tmp}" ] || ! grep -q "^>" "${out_tmp}"; then
        log_info "  → Downloading COMPLETE proteome (TrEMBL)"
        curl -s -G "https://rest.uniprot.org/uniprotkb/stream" \
            --data-urlencode "query=organism_id:${taxid}" \
            --data-urlencode "format=fasta" \
            --compressed --max-time $TIMEOUT \
            -o "${out_tmp}" 2>>"$ERROR_LOG"
        grep -q "^>" "${out_tmp}" 2>/dev/null && proteome_type="COMPLETE"
    fi

    if [ ! -s "${out_tmp}" ] || ! grep -q "^>" "${out_tmp}"; then
        rm -f "${out_tmp}"
        echo -e "${organism}\t${taxid}\t${proteome_id}\tNO_PROTEOME\t0\t-\tFAILED\tUniProt" >> "${PROTEOME_MANIFEST}"
        return 1
    fi

    mv "${out_tmp}" "${out_fa}"
    seq_count=$(grep -c "^>" "${out_fa}")
    log_success "  ✓ [UniProt] ${seq_count} sequences (${proteome_type}) → ${out_fa}"
    echo -e "${organism}\t${taxid}\t${proteome_id}\t${proteome_type}\t${seq_count}\t${out_fa}\tSUCCESS\tUniProt" >> "${PROTEOME_MANIFEST}"
    return 0
}

# --- NCBI fallback ---
download_proteome_ncbi() {
    local organism="$1"
    local taxid="$2"
    local out_fa="${PROTEOME_DIR}/${taxid}.fasta"
    local out_tmp="${out_fa}.tmp"

    log_info "  [NCBI] Downloading RefSeq proteins for TaxID ${taxid}..."

    esearch -db protein -query "txid${taxid}[Organism] AND refseq[filter]" | \
    efetch -format fasta > "${out_tmp}"

    if [ ! -s "${out_tmp}" ] || ! grep -q "^>" "${out_tmp}"; then
        log_error "  [NCBI] No RefSeq proteins found for TaxID ${taxid}"
        echo -e "${organism}\t${taxid}\t-\tNO_PROTEOME\t0\t-\tFAILED\tNCBI" >> "${PROTEOME_MANIFEST}"
        rm -f "${out_tmp}"
        return 1
    fi

    # Append TaxID and source tag to each header for downstream compatibility
    awk -v tid="${taxid}" \
        '/^>/ { print $0 " | taxid=" tid " | type=REFSEQ"; next } { print }' \
        "${out_tmp}" > "${out_fa}"
    rm -f "${out_tmp}"

    local seq_count
    seq_count=$(grep -c "^>" "${out_fa}")
    log_success "  ✓ [NCBI] ${seq_count} sequences (REFSEQ) → ${out_fa}"
    echo -e "${organism}\t${taxid}\trefseq_all\tREFSEQ\t${seq_count}\t${out_fa}\tSUCCESS\tNCBI" >> "${PROTEOME_MANIFEST}"
    return 0
}

# --- Per-organism orchestration ---
process_organism() {
    local organism="$1"
    organism=$(echo "$organism" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    [ -z "$organism" ] && return 0

    log_info "Processing: ${organism}"

    local taxid
    if ! taxid=$(get_taxid "$organism"); then
        return 1
    fi

    if download_proteome_uniprot "$organism" "$taxid"; then
        return 0
    fi

    log_info "  [Fallback] UniProt failed — trying NCBI..."
    if download_proteome_ncbi "$organism" "$taxid"; then
        return 0
    fi

    log_error "  Both UniProt and NCBI failed for: $organism"
    return 1
}

# --- Final statistics report ---
generate_statistics() {
    log_info "Generating statistics..."

    local total_orgs taxid_found taxid_notfound taxid_ambiguous
    local proteome_success proteome_failed
    local reference_count reviewed_count complete_count

    total_orgs=$(wc -l < "$ORGANISM_LIST")
    taxid_found=$(grep -c $'\tFOUND\t'     "$TAXID_MANIFEST" || echo 0)
    taxid_notfound=$(grep -c $'\tNOT_FOUND\t' "$TAXID_MANIFEST" || echo 0)
    taxid_ambiguous=$(grep -c $'\tAMBIGUOUS\t' "$TAXID_MANIFEST" || echo 0)
    proteome_success=$(grep -c $'\tSUCCESS\t'  "$PROTEOME_MANIFEST" || echo 0)
    proteome_failed=$(grep -c $'\tFAILED\t'   "$PROTEOME_MANIFEST" || echo 0)
    reference_count=$(grep -c $'\tREFERENCE\t' "$PROTEOME_MANIFEST" || echo 0)
    reviewed_count=$(grep -c $'\tREVIEWED\t'  "$PROTEOME_MANIFEST" || echo 0)
    complete_count=$(grep -c $'\tCOMPLETE\t'  "$PROTEOME_MANIFEST" || echo 0)

    {
        echo "========================================================================"
        echo "DOWNLOAD STATISTICS"
        echo "========================================================================"
        echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
        echo ""
        echo "ORGANISMS:"
        echo "  Total in list           : $total_orgs"
        echo ""
        echo "TAXID RESOLUTION (NCBI):"
        echo "  Resolved                : $taxid_found"
        echo "  Not found               : $taxid_notfound"
        echo "  Ambiguous               : $taxid_ambiguous"
        echo ""
        echo "PROTEOME DOWNLOADS:"
        echo "  Successful              : $proteome_success"
        echo "  Failed                  : $proteome_failed"
        echo ""
        echo "PROTEOME TYPES (UniProt only):"
        echo "  Reference               : $reference_count"
        echo "  Reviewed (Swiss-Prot)   : $reviewed_count"
        echo "  Complete (TrEMBL)       : $complete_count"
        echo ""
        echo "OUTPUT FILES:"
        echo "  FASTA files             : ${PROTEOME_DIR}/"
        echo "  TaxID manifest          : ${TAXID_MANIFEST}"
        echo "  Proteome manifest       : ${PROTEOME_MANIFEST}"
        echo "  Error log               : ${ERROR_LOG}"
        echo "========================================================================"
    } | tee "$STATS_FILE"
}

# --- Entry point ---
main() {
    local start_time end_time duration
    start_time=$(date +%s)

    init_environment
    check_dependencies

    if [ ! -f "${ORGANISM_LIST}" ]; then
        log_error "Organism list not found: ${ORGANISM_LIST}"
        exit 1
    fi

    log_info "Loading organism list into memory..."
    mapfile -t organisms_array < "${ORGANISM_LIST}"
    local total_lines=${#organisms_array[@]}
    log_info "Total organisms to process: ${total_lines}"
    echo "------------------------------------------------------------------------"

    local count=0 success=0 failed=0
    for organism in "${organisms_array[@]}"; do
        count=$((count + 1))
        echo ""
        log_info "Progress: ${count}/${total_lines}"
        [ -z "$organism" ] && continue

        if process_organism "$organism"; then
            success=$((success + 1))
        else
            failed=$((failed + 1))
        fi
        echo "------------------------------------------------------------------------"
    done

    end_time=$(date +%s)
    duration=$((end_time - start_time))

    echo ""
    echo "========================================================================"
    log_success "PROCESSING COMPLETE"
    echo "========================================================================"
    log_info "Elapsed time : ${duration} seconds"
    log_info "Successful   : ${success} | Failed: ${failed}"
    echo ""

    generate_statistics
}

main "$@"
