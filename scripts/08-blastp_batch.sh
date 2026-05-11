#!/usr/bin/env bash
# =============================================================================
# Script:      blastp_batch.sh
# Description: Runs BLASTp against the custom protein database for all query
#              FASTA files in the input directory. Applies post-search filters
#              for minimum percent identity and query coverage.
#
# Usage:       bash blastp_batch.sh
#              (edit the configuration block below if needed)
#
# Input:       input/blast/*.{fasta,fa,faa,txt}  — one file per gene family
#              blastDB                             — custom BLAST database
#
# Output (per query):
#   results/blast/output/<query>_blastp_results.tsv  — all hits (15-column tabular)
#   results/blast/filtered/<query>_blastp_filtered.tsv — hits passing identity/coverage
#   results/blast/log/<query>_blastp.log          — stderr log (if non-empty)
#
# Output columns: qseqid sseqid pident length mismatch gapopen qstart qend
#                 sstart send evalue bitscore qcovs qlen slen
#
# Dependencies: blastp (BLAST+ 2.17.0)
#
# Parameter rationale:
#   BLOSUM45   — optimized for cross-kingdom sequences with <40% pairwise
#                identity (inter-kingdom comparisons: fungi, bacteria, plants)
#   word_size 3 — increases sensitivity for divergent homologs
#   e-value 1e-5 — stringent statistical threshold; combined with post-search
#                  identity/coverage filters for two-stage quality control
#   SEG filter — masks low-complexity regions to suppress spurious alignments
# =============================================================================

# --- Configuration ---
INPUT_DIR="input/blast"
OUTPUT_DIR="results/blast"
OUTPUT_SUBDIR_BLASTP="${OUTPUT_DIR}/output"
OUTPUT_SUBDIR_FILTERED="${OUTPUT_DIR}/filtered"
OUTPUT_SUBDIR_LOG="${OUTPUT_DIR}/log"
DATABASE="data/blastDB"

EVALUE="0.00001"
WORD_SIZE="3"
MATRIX="BLOSUM45"
SEG="yes"

MIN_IDENTITY="30"
MIN_QCOV="50"

OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen"

# --- Pre-flight checks ---
echo "============================================"
echo "BLASTp batch search"
echo "============================================"
echo ""
echo "Configuration:"
echo "  E-value        : $EVALUE"
echo "  Word size      : $WORD_SIZE"
echo "  Matrix         : $MATRIX"
echo "  SEG filter     : $SEG"
echo "  Min. identity  : ${MIN_IDENTITY}%"
echo "  Min. qcovs     : ${MIN_QCOV}%"
echo ""

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"; exit 1
fi

mkdir -p "$OUTPUT_DIR" "$OUTPUT_SUBDIR_BLASTP" "$OUTPUT_SUBDIR_FILTERED" "$OUTPUT_SUBDIR_LOG"

if [ ! -f "${DATABASE}.phr" ]; then
    echo "ERROR: BLAST database not found: $DATABASE"
    echo "  Expected files: .phr .pin .psq"; exit 1
fi

if ! command -v blastp &> /dev/null; then
    echo "ERROR: blastp not found in PATH"
    echo "  Install with: conda install -c bioconda blast=2.17.0"; exit 1
fi

# --- Batch processing ---
contador=0; exitosos=0; fallidos=0

for input_file in "$INPUT_DIR"/*.{fasta,fa,faa,txt}; do
    [ -f "$input_file" ] || continue
    ((contador++))

    basename_file=$(basename "$input_file")
    filename="${basename_file%.*}"

    output_file="${OUTPUT_SUBDIR_BLASTP}/${filename}_blastp_results.tsv"
    filtered_file="${OUTPUT_SUBDIR_FILTERED}/${filename}_blastp_filtered.tsv"
    log_file="${OUTPUT_SUBDIR_LOG}/${filename}_blastp.log"

    echo "----------------------------------------"
    echo "Query [$contador]: $basename_file"

    # Run BLASTp
    blastp \
        -query  "$input_file" \
        -db     "$DATABASE" \
        -out    "$output_file" \
        -evalue "$EVALUE" \
        -outfmt "$OUTFMT" \
        -word_size "$WORD_SIZE" \
        -matrix    "$MATRIX" \
        -seg       "$SEG" \
        2> "$log_file"

    if [ $? -eq 0 ]; then
        ((exitosos++))
        num_hits=$(wc -l < "$output_file")
        echo "  ✓ BLASTp completed — total hits: $num_hits"

        if [ $num_hits -gt 0 ]; then
            # Post-search filter: retain hits meeting identity and coverage thresholds
            echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs\tqlen\tslen" \
                > "$filtered_file"
            awk -v min_id="$MIN_IDENTITY" -v min_cov="$MIN_QCOV" \
                '$3 >= min_id && $13 >= min_cov' "$output_file" >> "$filtered_file"

            num_filtered=$(( $(wc -l < "$filtered_file") - 1 ))
            echo "  Filtered hits (id≥${MIN_IDENTITY}%, cov≥${MIN_QCOV}%): $num_filtered"

            if [ $num_filtered -eq 0 ]; then
                rm "$filtered_file"
                echo "  ⚠ No hits passed the filtering thresholds"
            fi
        else
            echo "  ⚠ No hits found for this query"
        fi

        # Remove empty log files
        [ -s "$log_file" ] && echo "  See log: $log_file" || rm "$log_file"

    else
        echo "  ✗ ERROR: BLASTp failed — see $log_file"
        ((fallidos++))
    fi
    echo ""
done

# --- Summary ---
echo "============================================"
echo "Run complete"
echo "============================================"
echo "Files processed : $contador"
echo "Successful      : $exitosos"
echo "Failed          : $fallidos"
echo ""
echo "Results saved in: $OUTPUT_DIR"
echo "  *_blastp_results.tsv  — all hits"
echo "  *_blastp_filtered.tsv — filtered hits"
echo "  *_blastp.log          — run logs (if warnings/errors)"
echo "============================================"

if [ $contador -eq 0 ]; then
    echo "WARNING: No input files found in $INPUT_DIR"
    echo "  Expected extensions: .fasta .fa .faa .txt"
    exit 1
fi

[ $fallidos -gt 0 ] && exit 1
exit 0
