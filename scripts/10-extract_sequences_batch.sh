#!/usr/bin/env bash
# =============================================================================
# Script:      extract_sequences_batch.sh
# Description: Extracts FASTA sequences from the custom database using subject
#              sequence IDs from filtered BLASTp result tables. Processes all
#              gene families (subdirectories) found under the input directory.
#              Also generates a consolidated FASTA per gene family containing
#              all unique sequences across queries.
#
# Usage:       bash extract_sequences_batch.sh
#              (edit the configuration block below if needed)
#
# Input:       results/blast/unique/<protein_name>  — filtered BLASTp tables
#              data/hormoneDB.fasta      — merged protein database
#              scripts/GetSeqsFromFasta.py           — sequence extraction helper
#
# Output:      results/blast/extracted/<family>/<query>_sequences.fasta — per-query FASTA
#              results/blast/extracted/<family>/<family>_ALL_sequences.fasta — consolidated
#
# Dependencies: python3, biopython (via GetSeqsFromFasta.py), sort
# =============================================================================

# --- Configuration ---
INPUT_DIR="results/blast/unique/<protein_name>"
OUTPUT_DIR="results/blast/extracted/<protein_name>"
MULTIFASTA="data/hormoneDB.fasta"
PYTHON_SCRIPT="scripts/GetSeqsFromFasta.py"
TEMP_DIR="temp/ids"

SUBJECT_COLUMN=1        # Column in TSV containing subject sequence IDs (sseqid)
REMOVE_DUPLICATES="yes" # Remove duplicate IDs before extraction
CREATE_CONSOLIDATED="yes" # Create a merged FASTA per gene family

# --- Pre-flight checks ---
echo "============================================"
echo "Batch sequence extraction"
echo "============================================"
echo ""

for check in "$INPUT_DIR" "$MULTIFASTA" "$PYTHON_SCRIPT"; do
    if [ ! -e "$check" ]; then
        echo "ERROR: Not found: $check"; exit 1
    fi
done

if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found in PATH"; exit 1
fi

if ! python3 -c "import Bio" 2>/dev/null; then
    echo "ERROR: BioPython not installed"
    echo "  Install with: pip install biopython"; exit 1
fi

mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

echo "Configuration:"
echo "  Input directory  : $INPUT_DIR"
echo "  Output directory : $OUTPUT_DIR"
echo "  Database FASTA   : $MULTIFASTA"
echo "  Remove duplicates: $REMOVE_DUPLICATES"
echo "  Create consolidated: $CREATE_CONSOLIDATED"
echo ""

# --- Counters ---
total_families=0; total_files=0; exitosos=0; fallidos=0; sin_hits=0

# --- Main loop: one iteration per gene family (subdirectory) ---
for hormone_dir in "$INPUT_DIR"/*/; do
    [ -d "$hormone_dir" ] || continue

    hormone_name=$(basename "$hormone_dir")
    ((total_families++))

    echo "========================================"
    echo "Gene family [$total_families]: $hormone_name"
    echo "========================================"

    hormone_output_dir="${OUTPUT_DIR}/${hormone_name}"
    mkdir -p "$hormone_output_dir"

    all_ids_file="${TEMP_DIR}/${hormone_name}_all_ids.txt"
    > "$all_ids_file"

    family_file_count=0

    # Process each filtered TSV for this gene family
    for tsv_file in "$hormone_dir"*.tsv; do
        [ -f "$tsv_file" ] || continue
        ((total_files++)); ((family_file_count++))

        filename=$(basename "${tsv_file%.tsv}")
        echo ""
        echo "  Query [$family_file_count]: $(basename $tsv_file)"

        # Extract subject IDs from the specified column, skip comment lines
        ids_file="${TEMP_DIR}/${hormone_name}_${filename}_ids.txt"
        awk -v col="$SUBJECT_COLUMN" '!/^#/ { print $col }' "$tsv_file" > "$ids_file"

        num_ids=$(wc -l < "$ids_file")
        if [ "$num_ids" -eq 0 ]; then
            echo "    ⚠ No IDs found (empty file or no hits)"
            ((sin_hits++)); rm "$ids_file"; continue
        fi
        echo "    IDs extracted: $num_ids"

        # Remove duplicates if requested
        if [ "$REMOVE_DUPLICATES" = "yes" ]; then
            ids_unique="${ids_file%.txt}_unique.txt"
            sort -u "$ids_file" > "$ids_unique"
            num_unique=$(wc -l < "$ids_unique")
            removed=$((num_ids - num_unique))
            [ $removed -gt 0 ] && echo "    Unique IDs: $num_unique (removed $removed duplicates)"
            mv "$ids_unique" "$ids_file"
        fi

        cat "$ids_file" >> "$all_ids_file"

        # Extract sequences using the Python helper
        output_fasta="${hormone_output_dir}/${filename}_sequences.fasta"
        error_log="${TEMP_DIR}/error_${filename}.log"

        python3 "$PYTHON_SCRIPT" \
            --multifasta "$MULTIFASTA" \
            --ids        "$ids_file" \
            --output     "$output_fasta" 2> "$error_log"

        if [ -f "$output_fasta" ]; then
            num_seqs=$(grep -c "^>" "$output_fasta" 2>/dev/null || echo 0)
            if [ "$num_seqs" -gt 0 ]; then
                echo "    ✓ Sequences extracted: $num_seqs → $output_fasta"
                ((exitosos++)); rm -f "$error_log"
            else
                echo "    ⚠ No matching sequences found in the database"
                ((sin_hits++))
                [ -s "$error_log" ] && head -5 "$error_log"
            fi
        else
            echo "    ✗ ERROR: Extraction failed — output file not created"
            ((fallidos++))
            [ -s "$error_log" ] && cat "$error_log"
        fi
    done

    # Build consolidated FASTA for the gene family
    if [ "$CREATE_CONSOLIDATED" = "yes" ] && [ -s "$all_ids_file" ]; then
        echo ""
        echo "  Building consolidated FASTA for $hormone_name..."

        all_ids_unique="${TEMP_DIR}/${hormone_name}_all_ids_unique.txt"
        sort -u "$all_ids_file" > "$all_ids_unique"

        num_total=$(wc -l < "$all_ids_file")
        num_unique=$(wc -l < "$all_ids_unique")
        echo "    Total IDs: $num_total  |  Unique: $num_unique"

        consolidated="${hormone_output_dir}/${hormone_name}_ALL_sequences.fasta"
        python3 "$PYTHON_SCRIPT" \
            --multifasta "$MULTIFASTA" \
            --ids        "$all_ids_unique" \
            --output     "$consolidated" 2>&1 | grep -v "^INFO:"

        if [ -f "$consolidated" ]; then
            num_cons=$(grep -c "^>" "$consolidated")
            echo "    ✓ Consolidated FASTA: $consolidated ($num_cons sequences)"
        else
            echo "    ✗ ERROR: Could not create consolidated FASTA"
        fi
    fi

    echo ""
    echo "  Gene family summary — files processed: $family_file_count"
    echo ""
done

# --- Final summary ---
echo "============================================"
echo "Run complete"
echo "============================================"
echo "Gene families processed : $total_families"
echo "TSV files processed     : $total_files"
echo "Successful extractions  : $exitosos"
echo "Files with no hits      : $sin_hits"
echo "Failed extractions      : $fallidos"
echo ""
echo "Results saved in: $OUTPUT_DIR"
echo "Temporary ID files in  : $TEMP_DIR"
echo "============================================"
echo ""

read -r -p "Remove temporary ID files? (y/n): " clean_temp
if [[ "$clean_temp" =~ ^[Yy]$ ]]; then
    rm -rf "$TEMP_DIR"
    echo "✓ Temporary files removed"
else
    echo "Temporary files kept in: $TEMP_DIR"
fi

[ $fallidos -gt 0 ] && exit 1
exit 0
