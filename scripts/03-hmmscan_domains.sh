#!/usr/bin/env bash
# =============================================================================
# Script:      hmmscan_domains.sh
# Description: Scans seed protein sequences against the full Pfam-A database
#              to discover the domain architecture of each protein. 
#
#              Uses --cut_tc (Trusted Cutoff), the most stringent Pfam threshold.
#              Only domains whose score exceeds the lowest score in the manually
#              curated Pfam seed alignment are reported, minimizing false
#              positive domain assignments.
#
# Usage:       bash hmmscan_domains.sh <Pfam-A.hmm> <seeds.fasta> <output_dir>
#
#   <Pfam-A.hmm>   — path to local Pfam-A HMM database (pre-indexed with hmmpress)
#   <seeds.fasta>  — multi-FASTA with all seed protein sequences
#   <output_dir>   — directory where results will be written
#
# Output:      <output_dir>/<protein>/
#                <protein>_tc.tbl              — per-sequence hit table (--tblout)
#                <protein>_tc.dom              — per-domain hit table (--domtblout)
#                <protein>_tc.log              — full hmmscan text output
#                <protein>_domains_summary.tsv — curated domain list
#                                               (protein, domain, accession,
#                                                E-value, score, bias,
#                                                hmm_from, hmm_to, seq_from, seq_to)
#              <output_dir>/ALL_proteins_domain_summary.tsv
#                                             — consolidated table, all proteins
#
# Threshold options (edit THRESHOLD below if needed):
#   --cut_tc   Trusted Cutoff    — most stringent  [used in the published study]
#   --cut_ga   Gathering threshold — recommended by Pfam, slightly less strict
#   --cut_nc   Noise Cutoff      — permissive, for exploratory scans only
#
# Prerequisites:
#   Pfam-A.hmm must be indexed before first use:
#     hmmpress <Pfam-A.hmm>
#
# Dependencies: hmmscan, hmmpress (HMMER 3.x), awk
# =============================================================================

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <Pfam-A.hmm> <seeds.fasta> <output_dir>"
    echo "Example: $0 hmmer_db/Pfam-A.hmm input/seeds/seed_proteins.fasta output/hmmscan_domains"
    exit 1
fi

PFAM_DB="$1"
INPUT_FASTA="$2"
OUTPUT_DIR="$3"

THRESHOLD="--cut_tc"   # Change to --cut_ga or --cut_nc if needed
THREADS=4              # CPU threads for hmmscan

# --- Pre-flight checks ---
if ! command -v hmmscan &>/dev/null; then
    echo "Error: hmmscan not found. Ensure HMMER 3.x is installed and in PATH."
    exit 1
fi

if [ ! -f "$PFAM_DB" ]; then
    echo "Error: Pfam-A database not found: '$PFAM_DB'"
    exit 1
fi

if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input FASTA not found: '$INPUT_FASTA'"
    exit 1
fi

# Index Pfam-A if not yet done
if [ ! -f "${PFAM_DB}.h3i" ]; then
    echo "Pfam-A index not found. Running hmmpress..."
    hmmpress "$PFAM_DB"
fi

mkdir -p "$OUTPUT_DIR"

N_SEQS=$(grep -c "^>" "$INPUT_FASTA")

echo "============================================"
echo "Domain discovery — hmmscan vs Pfam-A"
echo "============================================"
echo "  Database  : $PFAM_DB"
echo "  Input     : $INPUT_FASTA ($N_SEQS sequences)"
echo "  Output    : $OUTPUT_DIR"
echo "  Threshold : $THRESHOLD"
echo "  Threads   : $THREADS"
echo ""

# --- Consolidated output header ---
CONSOLIDATED="${OUTPUT_DIR}/ALL_proteins_domain_summary.tsv"
printf "protein\tdomain_name\taccession\tE-value\tscore\tbias\thmm_from\thmm_to\tseq_from\tseq_to\n" \
    > "$CONSOLIDATED"

# --- Split input FASTA into one file per protein ---
TEMP_DIR="${OUTPUT_DIR}/.tmp_seqs"
mkdir -p "$TEMP_DIR"

awk -v tmp="$TEMP_DIR" '
/^>/ {
    if (outfile) close(outfile)
    header = substr($0, 2)
    split(header, parts, " ")
    name = parts[1]
    gsub(/[^A-Za-z0-9_.-]/, "_", name)
    outfile = tmp "/" name ".fasta"
    print > outfile
    next
}
{ if (outfile) print > outfile }
' "$INPUT_FASTA"

# --- Main loop: one hmmscan per protein ---
SUCCESS_COUNT=0
NO_DOMAIN_COUNT=0
ERROR_COUNT=0
CURRENT=0

for SEQ_FILE in "$TEMP_DIR"/*.fasta; do
    [ -f "$SEQ_FILE" ] || continue
    ((CURRENT++))

    PROTEIN=$(basename "$SEQ_FILE" .fasta)
    PROTEIN_DIR="${OUTPUT_DIR}/${PROTEIN}"
    mkdir -p "$PROTEIN_DIR"

    TBL_OUT="${PROTEIN_DIR}/${PROTEIN}_tc.tbl"
    DOM_OUT="${PROTEIN_DIR}/${PROTEIN}_tc.dom"
    LOG_OUT="${PROTEIN_DIR}/${PROTEIN}_tc.log"
    SUMMARY="${PROTEIN_DIR}/${PROTEIN}_domains_summary.tsv"

    echo "[$CURRENT/$N_SEQS] Scanning: $PROTEIN"

    hmmscan \
        $THRESHOLD \
        --tblout    "$TBL_OUT" \
        --domtblout "$DOM_OUT" \
        --cpu       "$THREADS" \
        "$PFAM_DB"  "$SEQ_FILE" \
        > "$LOG_OUT" 2>&1

    if [ $? -ne 0 ]; then
        echo "  ❌ Error: hmmscan failed — see $LOG_OUT"
        ((ERROR_COUNT++))
        continue
    fi

    DOMAIN_COUNT=$(grep -vc "^#" "$DOM_OUT" 2>/dev/null || echo 0)

    if [ "$DOMAIN_COUNT" -eq 0 ]; then
        echo "  ⚠  No domains found at $THRESHOLD threshold"
        ((NO_DOMAIN_COUNT++))
        printf "protein\tdomain_name\taccession\tE-value\tscore\tbias\thmm_from\thmm_to\tseq_from\tseq_to\n" \
            > "$SUMMARY"
        continue
    fi

    # Write per-protein summary TSV
    printf "protein\tdomain_name\taccession\tE-value\tscore\tbias\thmm_from\thmm_to\tseq_from\tseq_to\n" \
        > "$SUMMARY"

    grep -v "^#" "$DOM_OUT" | awk -v prot="$PROTEIN" 'NF > 20 {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            prot, $1, $2, $7, $8, $9, $16, $17, $20, $21
    }' | tee -a "$CONSOLIDATED" >> "$SUMMARY"

    echo "  ✅ Domains found: $DOMAIN_COUNT"
    grep -v "^#" "$DOM_OUT" | awk 'NF > 20 {
        printf "     %-30s acc=%-12s E=%-10s pos=%s-%s\n", $1, $2, $7, $20, $21
    }' | sort -u

    ((SUCCESS_COUNT++))
    echo ""
done

rm -rf "$TEMP_DIR"

# --- Summary ---
echo "============================================"
echo "Run complete"
echo "============================================"
echo "  Proteins scanned    : $N_SEQS"
echo "  ✅ With domains     : $SUCCESS_COUNT"
echo "  ⚠  No domains found : $NO_DOMAIN_COUNT"
echo "  ❌ Errors           : $ERROR_COUNT"
echo ""
echo "  Per-protein results : $OUTPUT_DIR/<protein>/"
echo "  Consolidated table  : $CONSOLIDATED"
echo ""
echo "Next step: use domain accessions from $CONSOLIDATED"
echo "as input list for extract_hmm_ids.sh"

[ $ERROR_COUNT -gt 0 ] && exit 1
exit 0