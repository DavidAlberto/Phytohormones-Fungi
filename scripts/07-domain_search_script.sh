#!/usr/bin/env bash
# =============================================================================
# Script:      domain_search_script.sh
# Description: Full domain search and extraction pipeline for a single Pfam
#              domain. Verifies that candidate sequences contain the expected
#              domain, extracts domain regions, aligns them to the HMM profile,
#              and optionally generates sequence logos.
#
#              This script was used to validate that candidate homologs possess
#              the complete domain architecture of the reference protein —
#              sequences lacking any expected domain were excluded from
#              downstream phylogenetic analysis.
#
# Usage:       bash domain_search_script.sh <PFAM_ID> <sequences.fasta> <output_dir>
#
#   <PFAM_ID>          — Pfam accession (e.g., PF00117)
#   <sequences.fasta>  — candidate protein sequences to evaluate
#   <output_dir>       — directory for output files
#
# Output (in output_<PFAM_ID>/):
#   <PFAM_ID>_hits.dom           — domain hit table (--domtblout format)
#   domains_<PFAM_ID>.fasta      — extracted domain sequences
#   coordinates_<PFAM_ID>.txt   — domain coordinates per sequence
#   <PFAM_ID>_alignment.fasta   — domain sequences aligned to the HMM
#   <PFAM_ID>_logo_*.png        — sequence logos (if WebLogo/HMMlogo available)
#
# Dependencies: hmmsearch, hmmfetch, hmmalign, esl-sfetch, esl-reformat
#              (HMMER 3.x / Easel); weblogo and hmmlogo are optional.
# =============================================================================

if [ $# -ne 3 ]; then
    echo "Usage: $0 <PFAM_ID> <sequences.fasta> <output_dir>"
    echo "Example: $0 PF00071 my_sequences.fasta results/hmmer/PF00071"
    exit 1
fi

PFAM_ID="$1"
SEQUENCES="$2"
OUTPUT_DIR="$3"

if [ ! -f "$SEQUENCES" ]; then
    echo "Error: Sequence file not found: $SEQUENCES"
    exit 1
fi

echo "=== Domain search pipeline for $PFAM_ID ==="
echo "Input : $SEQUENCES"
echo "Output: $OUTPUT_DIR"
echo ""

mkdir -p "$OUTPUT_DIR"

# --- Step 1: Fetch HMM profile from local Pfam-A database ---
echo "Step 1: Fetching HMM profile for $PFAM_ID..."
HMM_FILE="${OUTPUT_DIR}/${PFAM_ID}.hmm"

hmmfetch hmmer_db/Pfam-A.hmm "$PFAM_ID" > "$HMM_FILE"
if [ $? -ne 0 ]; then
    echo "  Error: Could not extract $PFAM_ID from Pfam-A.hmm"
    exit 1
fi
echo "  ✓ HMM profile saved: $HMM_FILE"

# --- Step 2: Search for the domain in candidate sequences ---
echo ""
echo "Step 2: Searching domain $PFAM_ID in sequences..."
HITS_FILE="${OUTPUT_DIR}/${PFAM_ID}_hits.dom"

hmmsearch --domtblout "$HITS_FILE" "$HMM_FILE" "$SEQUENCES" \
    > "${OUTPUT_DIR}/${PFAM_ID}_search.log"
if [ $? -ne 0 ]; then
    echo "  Error: hmmsearch failed"
    exit 1
fi

HITS_COUNT=$(grep -v "^#" "$HITS_FILE" | wc -l)
echo "  ✓ Search complete. Domain hits found: $HITS_COUNT"

if [ "$HITS_COUNT" -eq 0 ]; then
    echo "  ⚠ No $PFAM_ID domains found in the input sequences"
    exit 0
fi

# --- Step 3: Hit statistics ---
echo ""
echo "Step 3: Hit statistics..."
echo "  Sequences analyzed          : $(grep -c '^>' "$SEQUENCES")"
echo "  Domain hits                 : $HITS_COUNT"
echo "  Unique sequences with domain: $(grep -v '^#' "$HITS_FILE" | cut -f1 | sort -u | wc -l)"
echo "  Top 5 hits (by score):"
grep -v "^#" "$HITS_FILE" | sort -k13,13nr | head -5 | while read -r line; do
    target=$(echo "$line" | cut -d' ' -f1)
    score=$(echo  "$line" | cut -d' ' -f13)
    evalue=$(echo "$line" | cut -d' ' -f12)
    start=$(echo  "$line" | cut -d' ' -f20)
    end=$(echo    "$line" | cut -d' ' -f21)
    echo "    $target  score=$score  E-value=$evalue  pos=${start}-${end}"
done

# --- Step 4: Extract domain coordinates and sequences ---
echo ""
echo "Step 4: Extracting domain regions..."
COORDS_FILE="${OUTPUT_DIR}/coordinates_${PFAM_ID}.txt"
DOMAINS_FASTA="${OUTPUT_DIR}/domains_${PFAM_ID}.fasta"

awk '$1 !~ /^#/ && $1 != "" { print $1, $20, $21 }' "$HITS_FILE" > "$COORDS_FILE"
> "$DOMAINS_FASTA"

EXTRACTED=0
while read -r seq start end; do
    if [[ $start =~ ^[0-9]+$ ]] && [[ $end =~ ^[0-9]+$ ]] && [ "$start" -le "$end" ]; then
        esl-sfetch --range "${start}..${end}" "$SEQUENCES" "$seq" \
            >> "$DOMAINS_FASTA" 2>/dev/null && ((EXTRACTED++))
    fi
done < "$COORDS_FILE"

echo "  ✓ Domain sequences extracted: $EXTRACTED"

if [ ! -s "$DOMAINS_FASTA" ]; then
    echo "  ⚠ No domain sequences could be extracted"
    exit 1
fi

# --- Step 5: Optional sequence logos ---
echo ""
echo "Step 5: Generating sequence logos (optional tools)..."

LOGO_SEQS="${OUTPUT_DIR}/${PFAM_ID}_logo_from_sequences.png"
if command -v weblogo &>/dev/null; then
    weblogo -f "$DOMAINS_FASTA" -o "$LOGO_SEQS" -F png \
            --composition equiprobable \
            --title "Domain $PFAM_ID — logo from sequences" &>/dev/null \
        && echo "  ✓ Sequence logo: $LOGO_SEQS" \
        || echo "  ⚠ WebLogo failed"
else
    echo "  ⚠ weblogo not found — skipping sequence logo"
fi

LOGO_HMM="${OUTPUT_DIR}/${PFAM_ID}_logo_from_hmm.png"
if command -v hmmlogo &>/dev/null; then
    hmmlogo -H "$HMM_FILE" -o "$LOGO_HMM" &>/dev/null \
        && echo "  ✓ HMM logo: $LOGO_HMM" \
        || echo "  ⚠ hmmlogo failed"
else
    echo "  ⚠ hmmlogo not found — skipping HMM logo"
fi

# --- Step 6: Multiple sequence alignment against the HMM profile ---
echo ""
echo "Step 6: Aligning domain sequences to HMM profile..."
ALIGNMENT_STO="${OUTPUT_DIR}/${PFAM_ID}_alignment.sto"
ALIGNMENT_FASTA="${OUTPUT_DIR}/${PFAM_ID}_alignment.fasta"

hmmalign "$HMM_FILE" "$DOMAINS_FASTA" > "$ALIGNMENT_STO" 2>/dev/null
if [ $? -eq 0 ]; then
    esl-reformat fasta "$ALIGNMENT_STO" > "$ALIGNMENT_FASTA"
    echo "  ✓ Alignment (Stockholm) : $ALIGNMENT_STO"
    echo "  ✓ Alignment (FASTA)     : $ALIGNMENT_FASTA"
else
    echo "  ⚠ hmmalign failed"
fi

# --- Summary ---
echo ""
echo "=== SUMMARY ==="
echo "Domain        : $PFAM_ID"
echo "Hits found    : $HITS_COUNT"
echo "Output files in ${OUTPUT_DIR}/:"
echo "  Hit table          : $(basename "$HITS_FILE")"
echo "  Domain sequences   : $(basename "$DOMAINS_FASTA")"
echo "  Coordinates        : $(basename "$COORDS_FILE")"
[ -f "$ALIGNMENT_FASTA" ] && echo "  Alignment (FASTA)  : $(basename "$ALIGNMENT_FASTA")"
[ -f "$LOGO_SEQS" ]       && echo "  Logo (sequences)   : $(basename "$LOGO_SEQS")"
[ -f "$LOGO_HMM"  ]       && echo "  Logo (HMM)         : $(basename "$LOGO_HMM")"
echo ""
echo "✓ Pipeline completed successfully."
