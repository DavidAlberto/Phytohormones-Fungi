#!/usr/bin/env bash
# =============================================================================
# Script:      phylo_pipeline.sh
# Description: End-to-end phylogenetic analysis pipeline. Processes all FASTA
#              files in the input directory through the following steps:
#
#                1. Header renaming  — Genus_Species_TaxID_UniprotID_GeneName
#                2. Deduplication    — seqkit rmdup (exact sequence duplicates)
#                3. Clustering       — CD-HIT at 99% identity
#                4. Alignment        — MAFFT --auto --reorder
#                5. Guide tree       — FastTree -lg -gamma (approximate ML)
#                6. Alignment refine — PRANK -protein -iterate=3 (phylogeny-aware)
#                7. Trimming         — trimAl -automated1
#                8. Phylogeny        — IQ-TREE3 MFP + LG+C60, UFBoot 1000,
#                                      SH-aLRT 1000, -bnni, --symtest
#
#              Each step writes to a numbered output directory and implements
#              a checkpoint: if the expected output already exists, that step
#              is skipped. This allows safe resumption after interruption.
#
# Usage:       bash phylo_pipeline.sh
#              (edit the USER CONFIGURATION section below before running)
#
# Input:       01_sequences/*.fasta  — protein sequences per gene family
# Output:      02_filtering/  03_mafft/  04_fasttree/  05_prank/
#              06_trimal/     07_iqtree/
#
# Dependencies: seqkit, cd-hit, mafft, FastTree, prank, trimal, iqtree3
#              All tools must be available in the activated Conda environment.
# =============================================================================

# =============================================================================
# USER CONFIGURATION — edit these variables before running
# =============================================================================

# Name of the Conda environment containing all pipeline dependencies.
readonly CONDA_ENV_NAME="hormone"

# Base directory of your Conda installation.
# Common paths: ~/miniconda3  ~/anaconda3  /opt/conda  /usr/local
readonly CONDA_BASE="$HOME/miniconda3"

# Directory containing the input FASTA files (one per gene family).
readonly INPUT_DIR="results/phylo/01_sequences"

# =============================================================================
# END OF USER CONFIGURATION — do not edit below this line
# =============================================================================

# --- Conda initialization and environment activation ---
__conda_setup="$("${CONDA_BASE}/bin/conda" 'shell.bash' 'hook' 2>/dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
        . "${CONDA_BASE}/etc/profile.d/conda.sh"
    else
        echo "ERROR: Conda not found at '${CONDA_BASE}'."
        echo "  Edit the CONDA_BASE variable in the USER CONFIGURATION section."
        exit 1
    fi
fi
unset __conda_setup

conda activate "${CONDA_ENV_NAME}"

# --- Strict error handling (after conda activation) ---
set -e
set -u
set -o pipefail

# --- Directory setup ---
echo "[$(date)] Setting up output directories..."
mkdir -p "$INPUT_DIR/clean"
mkdir -p 02_filtering 03_mafft 04_fasttree 05_prank 06_trimal 07_iqtree

# --- Pipeline ---
shopt -s nullglob
FOUND_FILES=0

for INPUT_FILE in "$INPUT_DIR"/*.fasta; do

    FILENAME=$(basename -- "$INPUT_FILE")
    BASENAME="${FILENAME%.*}"

    # Checkpoint: skip the entire sample if IQ-TREE output already exists
    IQTREE_DIR="07_iqtree/${BASENAME}"
    IQTREE_FINAL_FILE="${IQTREE_DIR}/${BASENAME}.treefile"

    if [ -f "$IQTREE_FINAL_FILE" ]; then
        echo "=================================================================="
        echo "✅ ${BASENAME}: IQ-TREE output exists. Skipping."
        echo "=================================================================="
        FOUND_FILES=1
        continue
    fi

    echo "=================================================================="
    echo "▶ Processing: ${BASENAME}  [${FILENAME}]"
    echo "=================================================================="
    FOUND_FILES=1

    # --- Step 1: Header renaming ---
    # Renames headers to Genus_Species_TaxID_UniprotID_GeneName format.
    # Extracted fields: TaxID (OX=), species (OS=), gene name (GN=),
    # UniProt ID (third pipe-delimited field). Falls back to the raw header
    # if extraction fails. Appends a numeric suffix for duplicate names.
    echo "[$(date)] Step 1: Renaming FASTA headers..."
    CLEAN_FASTA="$INPUT_DIR/clean/${BASENAME}_clean.fasta"

    awk '
    /^>/ {
        full_header = $0;
        taxid=""; species=""; genename=""; uniprot_id="";

        if (match(full_header, /OX[=_]([0-9]+)/, a))          taxid = a[1];
        else if (match(full_header, /taxid[=_]([0-9]+)/, b))  taxid = b[1];
        else if (match(full_header, /^>([0-9]+)\|/, c))       taxid = c[1];

        if (match(full_header, /OS[=_]([A-Za-z0-9]+)[_ ]([A-Za-z0-9]+)/, d))
            species = d[1] "_" d[2];
        else if (match(full_header, /\[([A-Za-z0-9]+) ([A-Za-z0-9]+)/, f))
            species = f[1] "_" f[2];

        if (match(full_header, /GN=([A-Za-z0-9_.\/-]+)/, e))
            genename = e[1];

        if (match(full_header, /^>[^|]+\|[^|]+\|([A-Z0-9]+)\|/, g))
            uniprot_id = g[1];

        if (species != "") {
            gsub(/ \(nom\. inval\.\)/, "", species);
            gsub(/[\[\]\(\)]/, "", species);
            gsub(/\//, "_", species);
        }
        if (genename != "") {
            gsub(/\//, "_", genename);
            gsub(/[\[\]\(\)]/, "", genename);
        }

        if (species != "" && taxid != "") {
            base_name = species "_" taxid;
            if (uniprot_id != "") base_name = base_name "_" uniprot_id;
            if (genename    != "") base_name = base_name "_" genename;
        } else {
            base_name = substr($0, 2);
        }

        gsub(" ",              "_",  base_name);
        gsub(/\//,             "_",  base_name);
        gsub(/[\[\]\(\):;,]/, "",   base_name);
        gsub(/__+/,            "_",  base_name);

        count[base_name]++;
        final_name = (count[base_name] > 1) ? base_name "_" count[base_name] : base_name;

        print ">" final_name;
        next;
    }
    { print }
    ' "$INPUT_FILE" > "$CLEAN_FASTA"

    echo "   → Clean FASTA: $CLEAN_FASTA"

    # --- Step 2: Filtering (deduplication + clustering at 99% identity) ---
    FILTERED_DIR="02_filtering/${BASENAME}"
    FINAL_FASTA="${FILTERED_DIR}/${BASENAME}_final.fasta"
    echo "[$(date)] Step 2: Filtering (seqkit rmdup + CD-HIT 99%)..."

    if [ -f "$FINAL_FASTA" ]; then
        echo "   → Output exists. Skipping."
    else
        mkdir -p "$FILTERED_DIR"
        SEQKIT_OUT="${FILTERED_DIR}/${BASENAME}_seqkit.fasta"
        DUP_FILE="${FILTERED_DIR}/${BASENAME}_duplicates.txt"

        seqkit rmdup -s -i "$CLEAN_FASTA" -o "$SEQKIT_OUT" -d "$DUP_FILE"
        cd-hit -i "$SEQKIT_OUT" -o "$FINAL_FASTA" -c 0.99 -n 5 -d 0 -T 4 -M 8000
        rm -f "${FINAL_FASTA}.clstr"

        echo "   → Input : $(grep -c "^>" "$CLEAN_FASTA") sequences"
        echo "   → Output: $(grep -c "^>" "$FINAL_FASTA") sequences"
    fi

    # --- Step 3: Multiple sequence alignment (MAFFT) ---
    MAFFT_OUT="03_mafft/${BASENAME}_mafft.fasta"
    echo "[$(date)] Step 3: MAFFT alignment..."

    if [ -f "$MAFFT_OUT" ]; then
        echo "   → Output exists. Skipping."
    else
        mafft --auto --reorder "$FINAL_FASTA" > "$MAFFT_OUT"
    fi

    # --- Step 4: Guide tree (FastTree) ---
    GUIDE_TREE="04_fasttree/${BASENAME}_guide.nwk"
    echo "[$(date)] Step 4: FastTree guide tree..."

    if [ -f "$GUIDE_TREE" ]; then
        echo "   → Output exists. Skipping."
    else
        FastTree -lg -gamma -nosupport "$MAFFT_OUT" > "$GUIDE_TREE"
    fi

    # --- Step 5: Phylogeny-aware alignment refinement (PRANK) ---
    PRANK_DIR="05_prank/${BASENAME}"
    PRANK_PREFIX="${PRANK_DIR}/${BASENAME}"
    PRANK_OUT="${PRANK_PREFIX}.best.fas"
    echo "[$(date)] Step 5: PRANK alignment refinement (this may take time)..."

    if [ -f "$PRANK_OUT" ]; then
        echo "   → Output exists. Skipping."
    else
        mkdir -p "$PRANK_DIR"
        prank -d="$MAFFT_OUT" \
              -t="$GUIDE_TREE" \
              -o="$PRANK_PREFIX" \
              -protein +F -termgap -iterate=3 -showtree -showanc -uselogs -shortnames
    fi

    # --- Step 6: Alignment trimming (trimAl) ---
    TRIMAL_OUT="06_trimal/${BASENAME}.trimal.fasta"
    TRIMAL_HTML="06_trimal/${BASENAME}_report.html"
    echo "[$(date)] Step 6: trimAl..."

    if [ -f "$TRIMAL_OUT" ]; then
        echo "   → Output exists. Skipping."
    else
        trimal -in "$PRANK_OUT" \
               -out "$TRIMAL_OUT" \
               -automated1 \
               -htmlout "$TRIMAL_HTML"
    fi

    # --- Step 7: Maximum-likelihood phylogenetic inference (IQ-TREE3) ---
    # Model selection: ModelFinder Plus (MFP) with additional profile mixture
    # models LG+C60 and LG+F+C60, which better capture compositional
    # heterogeneity across distantly related taxa.
    # Support: ultrafast bootstrap (UFBoot 1000) + SH-aLRT test (1000).
    IQTREE_PREFIX="${IQTREE_DIR}/${BASENAME}"
    echo "[$(date)] Step 7: IQ-TREE3 phylogenetic inference..."

    if [ -f "$IQTREE_FINAL_FILE" ]; then
        echo "   → Output exists. Skipping."
    else
        mkdir -p "$IQTREE_DIR"
        iqtree3 --prefix "$IQTREE_PREFIX" \
                -s "$TRIMAL_OUT" \
                -m MFP \
                -madd LG+C60,LG+F+C60 \
                -B 1000 \
                -alrt 1000 \
                -bnni \
                -T AUTO \
                --symtest

        # Copy treefile as .nwk for downstream visualization scripts
        cp "${IQTREE_PREFIX}.treefile" "${IQTREE_PREFIX}.nwk"
    fi

    echo "✅ ${BASENAME}: all steps completed."
done

shopt -u nullglob

if [ "$FOUND_FILES" -eq 0 ]; then
    echo "=================================================================="
    echo "WARNING: No .fasta files found in '${INPUT_DIR}'."
    echo "=================================================================="
    exit 1
fi

echo "================================================================================"
echo "[$(date)] ALL PHYLOGENETIC ANALYSES COMPLETED"
echo "================================================================================"
