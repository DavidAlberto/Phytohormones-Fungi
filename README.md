# When Fungi Speak the Language of Plants: Shared Phytohormones with Divergent Meanings

This repository describes a **reproducible bioinformatics workflow** designed to identify, curate, and analyze homologs of phytohormone-associated genes in filamentous fungi. The strategy integrates sequence similarity searches, phylogenetic reconstruction, and structural modeling to support evolutionary and functional inference.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/field-Bioinformatics-blue.svg)]()
[![Phylogenetics](https://img.shields.io/badge/analysis-Phylogenetics-green.svg)]()

---

## Table of Contents

* **[Overview](#overview)**
* **[Installation](#installation)**
   * [Conda Environment](#conda-environment)
   * [Docker Image](#docker-image)
* **[Data Requirements](#data-requirements)**
* **[Pipeline Steps](#pipeline-steps)**
   * [1. Construction of Custom Protein Database](#1-construction-of-a-custom-protein-database)
   * [2. Selection of Seed Sequences](#2-selection-of-seed-sequences)
   * [3. Homology Detection via BLASTp](#3-homology-detection-via-blastp)
   * [4. Filtering and Curation of BLAST Results](#4-filtering-and-curation-of-blast-results)
   * [5. Sequence Retrieval](#5-sequence-retrieval)
   * [6. Multiple Sequence Alignment](#6-multiple-sequence-alignment)
   * [7. Phylogeny-Aware Alignment Refinement](#7-phylogeny-aware-alignment-refinement)
   * [8. Phylogenetic Inference](#8-phylogenetic-inference)
   * [9. Phylogenetic Tree Visualization](#9-phylogenetic-tree-visualization)
   * [10. Structural Modeling](#10-structural-modeling-complementary-analysis)
* **[Example Usage](#example-usage)**
* **[Directory Structure](#directory-structure)**
* **[Software Requirements](#software-requirements)**
* **[Reproducibility Considerations](#reproducibility-considerations)**
* **[Applications](#applications)**
* **[Best Practices](#best-practices)**
* **[Contributing](#contributing)**
* **[Citation](#citation)**
* **[License](#license)**
* **[Acknowledgments](#acknowledgments)**
* **[Contact](#contact)**
* **[FAQ](#faq)**

---

## Overview

The pipeline was conceived to:

- Select reliable **seed sequences** based on curated functional evidence
- Identify homologs using **BLASTp searches** against custom fungal databases
- **Filter and curate** BLAST results to select best hits per organism
- **Extract sequences** systematically from multifasta databases
- Generate high-quality multiple sequence alignments
- Infer robust phylogenetic relationships
- **Visualize phylogenies** with integrated taxonomic metadata
- Support structural interpretation through **AlphaFold** models

The workflow is platform-agnostic and can be adapted to any gene family of interest.

---

## Installation

### Conda Environment

We provide a Conda environment file with all required dependencies for the phylogenetic analysis pipeline.

> **💡 Tip:** If you prefer not to manage dependencies, use our [Docker image](#docker-image) which has everything pre-installed and tested.

**Option 1: Create environment from YAML file**

```bash
# Clone the repository
git clone https://github.com/yourusername/phytohormone-phylogenomics.git
cd phytohormone-phylogenomics

# Create the conda environment
conda env create -f environment.yml

# Activate the environment
conda activate hormone-phylo
```

**Option 2: Manual environment creation**

Create the environment with the exact versions we tested:

```bash
# Create a new conda environment
conda create -n hormone-phylo python=3.12

# Activate the environment
conda activate hormone-phylo

# Install bioinformatics tools with tested versions
conda install -c bioconda -c conda-forge \
  blast=2.17.0 \
  mafft=7.526 \
  fasttree=2.1.11 \
  trimal=1.5.0 \
  iqtree=3.0.1 \
  seqkit=2.10.1 \
  cd-hit=4.8.1 \
  hmmer=3.4 \
  prank \
  biopython \
  numpy
```

**Verify installation:**

```bash
# Check versions of key tools
blastp -version          # Should show BLAST 2.17.0+
mafft --version          # Should show v7.526
iqtree3 --version        # Should show IQ-TREE multicore version 3.0.1
trimal --version         # Should show trimAl v1.5.0
```

---

### Docker Image

We provide a Docker image with all dependencies pre-installed for maximum reproducibility and portability across different systems. The image contains all software versions tested and validated for this pipeline.

**Pull the Docker image:**

```bash
docker pull davidalbertoge/hormone-analysis:latest
```

**Run the container:**

```bash
# Run interactively with current directory mounted
docker run -v $(pwd):/home/ -it davidalbertoge/hormone-analysis:latest
```

**What's included in the Docker image:**
- All bioinformatics tools with tested versions
- Pre-configured environment ready to run

**Build from Dockerfile (optional):**

```bash
# Build the image locally from provided Dockerfile
docker build -t hormone-analysis:local .

# Run your local build
docker run -v $(pwd):/home/ -it hormone-analysis:local
```

---

## Data Requirements

### Input Data

**1. Proteome Databases**

The pipeline requires protein sequence databases from target organisms. We recommend the following sources:

| Database | Source | Format | URL |
|----------|--------|--------|-----|
| **UniProt** | Reviewed proteins (Swiss-Prot) | FASTA | https://www.uniprot.org/downloads |
| **NCBI RefSeq** | Reference protein sequences | FASTA | https://www.ncbi.nlm.nih.gov/refseq/ |
| **JGI MycoCosm** | Fungal genomes | FASTA | https://mycocosm.jgi.doe.gov/ |
| **FungiDB** | Fungal genomics resource | FASTA | https://fungidb.org/ |

**2. Seed Sequences**

Experimentally characterized protein sequences used as BLAST queries:

- Should be in FASTA format
- One file per gene family or multiple families in separate files
- Headers should contain gene names and organism information

**Example seed sequence format:**

```
>sp|P12345|GENE_ORGANISM Gene description OS=Organism name OX=12345 GN=geneName
MAKTPVQIWSFLKDHGFSDKHGFKJHGFKJDHGFKJDHGF...
```

**3. Taxonomic Metadata**

For phylogenetic tree visualization, prepare a metadata file with:

| Column | Description | Example |
|--------|-------------|---------|
| TaxID | NCBI Taxonomy ID | 5476 |
| Organism | Scientific name | Candida albicans |
| Kingdom | Taxonomic kingdom | Fungi |
| Phylum | Taxonomic phylum | Ascomycota |
| EarlyDivergent | Basal lineage flag | TRUE/FALSE |

**Metadata file format (TSV):**

```tsv
TaxID	Organism	Kingdom	Phylum	EarlyDivergent
5476	Candida albicans	Fungi	Ascomycota	FALSE
367110	Neurospora crassa	Fungi	Ascomycota	FALSE
284811	Aspergillus fumigatus	Fungi	Ascomycota	FALSE
```

### AlphaFold Requirements

If using AlphaFold for structural modeling:

**AlphaFold web server**: https://alphafold.ebi.ac.uk/

### ChimeraX Requirements

For structural visualization:

- **Download**: https://www.cgl.ucsf.edu/chimerax/download.html
- **Platform**: Available for Linux, macOS, and Windows
- **License**: Free for non-commercial use
- **System**: OpenGL 3.3 or later required

---

## Pipeline Steps

### 1. Construction of a Custom Protein Database

Proteomes from selected organisms were integrated into a **non-redundant multifasta database** representing the taxonomic diversity relevant to the study. This database constituted the search space for downstream homology detection.

**Database construction:**

```bash
# Concatenate individual proteome files into a single multifasta
cat organism1.fasta organism2.fasta organism3.fasta > fungal_proteomes.fasta

# Create BLAST database
makeblastdb -in fungal_proteomes.fasta -dbtype prot -out fungal_db
```

**Key considerations:**
- Headers should follow a consistent format (e.g., `TaxID|source|UniProtID|...`)
- Document the taxonomic composition of the database

---

### 2. Selection of Seed Sequences

Initial queries were chosen based on:

- **Sequence homology** to experimentally characterized proteins
- Presence of expected functional features described in the literature
- Structural conservation of key catalytic or regulatory domains

This step ensured that only biologically meaningful references were used to initiate similarity searches, minimizing propagation of misannotated sequences.

**Criteria for seed selection:**
- Proteins with experimental validation (in literature)
- Complete sequences with well-defined domain architecture
- Representatives from diverse taxonomic groups when applicable

---

### 3. Homology Detection via BLASTp

Candidate homologs were identified using BLASTp against the custom database with parameters optimized for divergent fungal proteins.

**BLASTp execution:**

```bash
blastp \
  -query seeds.fasta \
  -db fungal_db \
  -out blastp_results.tsv \
  -evalue 0.00001 \
  -word_size 3 \
  -matrix BLOSUM45 \
  -seg yes \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen"
```

**Parameter rationale:**
- **E-value threshold (0.00001)**: Stringent cutoff for statistical significance
- **Word size (3)**: Increased sensitivity for detecting divergent homologs
- **Matrix (BLOSUM45)**: Optimized for evolutionarily distant sequences
- **SEG filter (yes)**: Masks low-complexity regions to reduce spurious hits

**Output format:**
- Extended tabular format (outfmt 6) with 15 columns
- Includes query coverage (`qcovs`) and sequence lengths (`qlen`, `slen`)

---

### 4. Filtering and Curation of BLAST Results

Raw BLAST results were systematically filtered to retain only high-quality hits and eliminate redundancy. The filtering process addressed a critical challenge: **selecting the best hit per organism** when multiple paralogs or isoforms were detected.

**Filtering strategy:**

```r
# Filtering criteria applied:
# 1. Minimum identity: ≥30%
# 2. Minimum query coverage: ≥50%
# 3. Organism extraction from subject IDs
# 4. Selection of best hit per Query-Organism pair

# Prioritization for best hit selection:
# 1. Bit score (highest) - alignment quality
# 2. Query coverage (highest) - % of query aligned
# 3. Percent identity (highest) - sequence similarity
# 4. E-value (lowest) - statistical significance
```

**Key innovation:**
- Organism IDs were extracted from subject sequence headers (TaxID or species code)
- Multiple hits from the same organism were consolidated to a **single best representative**
- This approach ensures one sequence per organism per query, facilitating downstream phylogenetic analysis

**Quality control metrics:**
- Distribution of identity percentages
- Query coverage statistics
- E-value ranges
- Number of organisms per gene family

**Outputs:**
- Consolidated table with best hits per Query-Organism combination
- Individual tables for each query gene
- Statistical summary of filtering results
- Visualization of quality metrics (identity vs. coverage)

---

### 5. Sequence Retrieval

Selected hits were systematically retrieved from the multifasta database based on their sequence IDs.

**Extraction process:**

```bash
# Extract subject IDs from filtered BLAST results
awk '{print $2}' filtered_results.tsv > sequence_ids.txt

# Remove duplicate IDs
sort -u sequence_ids.txt > unique_ids.txt

# Extract sequences from database
# Using BioPython or similar tools
seqkit grep -f unique_ids.txt fungal_proteomes.fasta > extracted_sequences.fasta
```

**Quality checks:**
- Verify all IDs were successfully retrieved
- Check for truncated or incomplete sequences
- Confirm sequence lengths match expected protein sizes
- Remove sequences with extensive gaps or ambiguous residues

---

### 6. Multiple Sequence Alignment

High-quality MSAs were generated for each protein family using progressive alignment strategies.

**Initial alignment:**

```bash
mafft --auto --reorder curated_sequences.fasta > alignment.fasta
```

**MAFFT parameters:**
- `--auto`: Automatically selects alignment strategy based on dataset size
- `--reorder`: Orders sequences by similarity for improved visualization

**Alignment considerations:**
- Conservation of catalytic and regulatory motifs
- Proper alignment of domain boundaries
- Detection of divergent or potentially mispredicted regions
- Identification of regions requiring manual curation

---

### 7. Phylogeny-Aware Alignment Refinement

To improve alignment quality, a guide tree was used to inform positional homology inference.

**Guide tree generation:**

```bash
FastTree -lg -gamma -nosupport alignment.fasta > guide.nwk
```

**Refinement with PRANK:**

```bash
prank -d=alignment.fasta \
      -t=guide.nwk \
      -protein +F -termgap -iterate=3 \
      -showtree -showanc -uselogs -shortnames
```

**PRANK advantages:**
- Distinguishes insertions from deletions based on phylogeny
- Preserves insertion events as positional homologs
- Iteratively refines alignment (3 iterations recommended)

**Trimming poorly aligned regions:**

```bash
trimal -in prank.best.fas \
       -out alignment_trimmed.fasta \
       -automated1 \
       -htmlout report.html
```

**trimAl parameters:**
- `-automated1`: Heuristic method balancing alignment quality and information retention
- HTML report provides visual assessment of trimmed positions

---

### 8. Phylogenetic Inference

Gene trees were reconstructed using maximum-likelihood with extended model exploration and statistical support estimation.

**IQ-TREE execution:**

```bash
iqtree3 --prefix analysis \
        -s alignment_trimmed.fasta \
        -m MFP \
        -madd LG+C60,LG+F+C60 \
        -B 1000 \
        -alrt 1000 \
        -bnni \
        -T AUTO \
        --symtest
```

**IQ-TREE parameters:**
- `-m MFP`: ModelFinder Plus for automatic model selection
- `-madd LG+C60,LG+F+C60`: Include profile mixture models for heterogeneous data
- `-B 1000`: Ultrafast bootstrap with 1000 replicates
- `-alrt 1000`: SH-aLRT test for additional support assessment
- `-bnni`: Optimize bootstrap trees with NNI
- `--symtest`: Test for symmetry in substitution patterns

**Statistical support:**
- UFBoot ≥95%: Strong support
- SH-aLRT ≥80%: Additional confidence measure
- Both metrics reported at internal nodes

**Tree interpretation:**
The resulting topologies were used to:
- Distinguish orthologous clades
- Detect lineage-specific expansions or losses
- Propose functional diversification events
- Identify potential horizontal gene transfer candidates

---

### 9. Phylogenetic Tree Visualization

Phylogenetic trees were visualized with integrated taxonomic metadata to enhance biological interpretation.

**Visualization approach:**

```r
# Tree visualization with ggtree in R
# Features:
# - Cladogram or phylogram layouts
# - Bootstrap support values displayed at nodes
# - Branch lengths optionally shown
# - Tips colored by taxonomic Kingdom
# - Tips shaped by Phylum
# - Highlighting of specific clades

# Metadata integration:
# - TaxID-based matching with tip labels
# - Multiple matching strategies for robustness
# - Color schemes optimized for clarity
# - Phylum abbreviations for space efficiency
```

**Key visualization elements:**

1. **Layout options:**
   - Cladogram (topology-focused)
   - Phylogram (branch length-based)

2. **Statistical support:**
   - UFBoot/SH-aLRT values at internal nodes
   - Configurable size and positioning
   - Color-coded thresholds (optional)

3. **Taxonomic metadata:**
   - **Kingdom**: Color-coded labels (e.g., Fungi, Viridiplantae, Metazoa)
   - **Phylum**: Shape-coded tip points
   - **Early divergence**: Special notation for basal lineages

4. **Clade highlighting:**
   - Visual emphasis on taxonomic groups of interest
   - Configurable transparency and colors

---

### 10. Structural Modeling

Representative sequences from major clades were modeled with **AlphaFold** to:

- Evaluate conservation of active sites
- Compare structural features among paralogs
- Support functional hypotheses derived from phylogeny
- Identify potential cryptic paralogs based on structural divergence

**Structural analysis in ChimeraX:**

- Conservation of catalytic triads or active site residues
- Domain organization and relative positioning
- Structural basis for functional diversification

---

## Example Usage

This section provides a complete walkthrough of the pipeline using example data.

### Quick Start

```bash
# 1. Activate the conda environment
conda activate hormone-phylo

# 2. Navigate to your project directory
cd phytohormone-phylogenomics

# 3. Create directory structure
mkdir -p data/{proteomes,seeds,blast,filtered,sequences,alignments,trees}

# 4. Build BLAST database
makeblastdb -in data/proteomes/fungal_db.fasta -dbtype prot -out data/proteomes/fungal_db

# 5. Run BLAST search
blastp \
  -query data/seeds/IPT_seeds.fasta \
  -db data/proteomes/fungal_db \
  -out data/blast/IPT_results.tsv \
  -evalue 0.00001 \
  -word_size 3 \
  -matrix BLOSUM45 \
  -seg yes \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen"

# 6. Filter BLAST results (using R)
Rscript scripts/filter_blast_results.R \
  --input data/blast/IPT_results.tsv \
  --output data/filtered/IPT_filtered.tsv \
  --min-identity 30 \
  --min-coverage 50

# 7. Extract sequences
seqkit grep -f data/filtered/IPT_ids.txt \
  data/proteomes/fungal_db.fasta \
  > data/sequences/IPT_sequences.fasta

# 8. Multiple sequence alignment
mafft --auto --reorder \
  data/sequences/IPT_sequences.fasta \
  > data/alignments/IPT_mafft.fasta

# 9. Generate guide tree
FastTree -lg -gamma -nosupport \
  data/alignments/IPT_mafft.fasta \
  > data/alignments/IPT_guide.nwk

# 10. Refine alignment with PRANK
prank -d=data/alignments/IPT_mafft.fasta \
      -t=data/alignments/IPT_guide.nwk \
      -o=data/alignments/IPT_prank \
      -protein +F -termgap -iterate=3

# 11. Trim alignment
trimal -in data/alignments/IPT_prank.best.fas \
       -out data/alignments/IPT_trimmed.fasta \
       -automated1

# 12. Phylogenetic inference
iqtree3 --prefix data/trees/IPT \
        -s data/alignments/IPT_trimmed.fasta \
        -m MFP \
        -madd LG+C60,LG+F+C60 \
        -B 1000 \
        -alrt 1000 \
        -bnni \
        -T AUTO

# 13. Visualize tree (using R)
Rscript scripts/visualize_tree.R \
  --tree data/trees/IPT.treefile \
  --metadata data/metadata.tsv \
  --output figures/IPT_tree.svg
```

### Example Output Structure

```
results/
├── IPT/
│   ├── 01_blast/
│   │   └── IPT_results.tsv
│   ├── 02_filtered/
│   │   ├── IPT_filtered.tsv
│   │   └── IPT_stats.txt
│   ├── 03_sequences/
│   │   └── IPT_sequences.fasta
│   ├── 04_alignments/
│   │   ├── IPT_mafft.fasta
│   │   ├── IPT_prank.best.fas
│   │   └── IPT_trimmed.fasta
│   ├── 05_trees/
│   │   ├── IPT.treefile
│   │   ├── IPT.log
│   │   └── IPT.iqtree
│   └── 06_figures/
│       ├── IPT_tree.svg
│       ├── IPT_tree.png
│       └── IPT_legend.svg
└── CKX/
    └── ... (similar structure)
```

---

## Directory Structure

Recommended project organization for reproducibility:

```
phytohormone-phylogenomics/
│
├── README.md
├── LICENSE
├── environment.yml          # Conda environment specification
├── Dockerfile               # Docker container specification
│
├── data/                    # Input data (not tracked in git)
│   ├── proteomes/           # Protein sequence databases
│   │   ├── organism1.fasta
│   │   ├── organism2.fasta
│   │   └── fungal_db.fasta  # Combined database
│   ├── seeds/               # Query sequences
│   │   ├── IPT_seeds.fasta
│   │   ├── CKX_seeds.fasta
│   │   └── ...
│   └── metadata.tsv         # Taxonomic metadata
│
├── scripts/                 # Analysis scripts
│   ├── filter_blast_results.R
│   ├── visualize_tree.R
│   ├── run_pipeline.sh
│   └── fetch_taxonomy.py
│
├── results/                 # Analysis outputs (not tracked in git)
│   ├── IPT/
│   ├── CKX/
│   └── ...
│
├── figures/                 # Publication-ready figures
│   ├── IPT_tree.svg
│   ├── IPT_tree.png
│   └── ...
│
└── docs/                    # Additional documentation
    ├── methods.md
    └── supplementary.md
```

**Git configuration (.gitignore):**

```
# Large data files
data/proteomes/
results/
*.fasta
*.tsv
*.nwk
*.treefile

# Temporary files
*.log
*.tmp
*.bak

# R files
.Rhistory
.RData
.Rproj.user/

# Python cache
__pycache__/
*.pyc
```

---

## Software Requirements

### Core Tools

The following table lists the **tested and validated versions** included in our Docker image and Conda environment:

| Tool | Version | Purpose | Installation |
|------|---------|---------|--------------|
| **BLAST+** | 2.17.0 | Homology searches | `conda install -c bioconda blast=2.17.0` |
| **MAFFT** | 7.526 | Multiple sequence alignment | `conda install -c bioconda mafft=7.526` |
| **FastTree** | 2.1.11 | Guide tree generation | `conda install -c bioconda fasttree` |
| **PRANK** | 170427 | Phylogeny-aware alignment refinement | `conda install -c bioconda prank` |
| **trimAl** | 1.5.0 | Alignment trimming | `conda install -c bioconda trimal=1.5.0` |
| **IQ-TREE** | 3.0.1 | Phylogenetic inference | `conda install -c bioconda iqtree=3.0.1` |
| **CD-HIT** | 4.8.1 | Sequence clustering | `conda install -c bioconda cd-hit` |
| **SeqKit** | 2.10.1 | Sequence manipulation | `conda install -c bioconda seqkit=2.10.1` |
| **HMMER** | 3.4 | Profile HMM searches (optional) | `conda install -c bioconda hmmer=3.4` |

**Optional tools:**

| Tool | Version | Purpose |
|------|---------|---------|
| AlphaFold | ≥2.3.0 | Protein structure prediction |
| ChimeraX | ≥1.5 | Structural visualization |

### R Packages

**Recommended R version:** ≥4.0.0

**Data manipulation and filtering:**
- `tidyverse` (≥1.3.0): Data wrangling and visualization
- `dplyr` (≥1.0.0): Data frame manipulation
- `stringr` (≥1.4.0): String operations

**Phylogenetic visualization:**
- `ggtree` (≥3.0.0): Tree visualization (Bioconductor)
- `treeio` (≥1.16.0): Tree I/O operations (Bioconductor)
- `ape` (≥5.5): Phylogenetic analysis
- `phangorn` (≥2.8.0): Additional phylogenetic tools

**Graphics and export:**
- `ggplot2` (≥3.3.0): Grammar of graphics
- `svglite` (≥2.1.0): SVG export
- `ggnewscale` (≥0.4.0): Multiple scales in ggplot2
- `here` (≥1.0.0): Project-relative paths

**Installation:**
```r
# CRAN packages
install.packages(c("tidyverse", "svglite", "here", "ape", "phangorn"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio"))
```

### System Requirements

**Recommended:**
- RAM: 16 GB (32 GB for large datasets)
- Disk space: 100 GB (databases can be large)
- CPU: 8+ cores (for parallel processing)
- SSD storage (significantly improves I/O performance)

---

## Reproducibility Considerations

The workflow relies exclusively on widely adopted tools in computational biology. We have tested and validated specific versions to ensure reproducibility across different systems.

### Version Control and Compatibility

- **Tested versions**: All tools have been validated with the versions listed in [Software Requirements](#software-requirements)
- **Docker container**: Pre-configured environment available at `davidalbertoge/hormone-analysis:latest`
- **Conda environment**: Reproducible environment specification provided in `environment.yml`
- **Parameter transparency**: All critical parameters explicitly stated in each pipeline step
- **Platform independence**: Runs on Linux, macOS, and Windows (with WSL2)

### Ensuring Reproducibility

1. **Use version-controlled software:**
   ```bash
   # Recommended: Use our Docker container
   docker pull davidalbertoge/hormone-analysis:latest
   
   # Or create Conda environment with pinned versions
   conda env create -f environment.yml
   ```

2. **Document your analysis:**
   - Record exact command lines used
   - Note any parameter changes from defaults
   - Save log files from each step
   - Version control your scripts with Git

3. **Random seed handling:**
   - Bootstrap analyses use random sampling
   - IQ-TREE uses `-seed` parameter for reproducibility
   - FastTree results may vary slightly between runs
   - Document random seeds when applicable

4. **Hardware considerations:**
   - Number of CPU threads can affect results (especially with stochastic methods)
   - Use `-T` parameter to control threads in IQ-TREE
   - Memory allocation can impact performance but not results

---

## Applications

This strategy is suitable for:

- **Evolutionary reconstruction** of phytohormone-related pathways in fungi
- **Identification of candidate functional homologs** for experimental validation
- **Comparative analysis** of fungal signaling networks across taxonomic scales
- **Detection of horizontal gene transfer** events
- **Functional prediction** based on phylogenetic context
- **Selection of targets** for structural or biochemical characterization

---

## Best Practices

### General Recommendations

1. **Seed sequence quality**: Use experimentally validated proteins when available
2. **Database curation**: Remove redundant sequences and verify annotations
3. **BLAST parameters**: Adjust word size and matrix based on evolutionary distance
4. **Alignment inspection**: Manually verify critical motifs and domains
5. **Model selection**: Allow IQ-TREE to explore multiple models
6. **Bootstrap support**: Report both UFBoot and SH-aLRT values
7. **Metadata integration**: Ensure taxonomic IDs are current and accurate

---

## Contributing

We welcome contributions to improve this pipeline! Here's how you can help:

### Reporting Issues

If you encounter bugs or have suggestions:

1. Check if the issue already exists in [Issues](https://github.com/DavidAlberto/Phytohormones-Fungi/issues)
2. Create a new issue with:
   - **Clear title** describing the problem
   - **Steps to reproduce** the issue
   - **Expected behavior** vs actual behavior
   - **Environment details** (OS, software versions)
   - **Error messages** or log files

### Submitting Changes

1. **Fork the repository**
   ```bash
   git clone https://github.com/yourusername/phytohormone-phylogenomics.git
   cd phytohormone-phylogenomics
   ```

2. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

3. **Make your changes**
   - Follow existing code style
   - Add comments explaining complex steps
   - Update documentation if needed

4. **Test your changes**
   ```bash
   # Run on test dataset
   bash scripts/test_pipeline.sh
   ```

5. **Commit and push**
   ```bash
   git add .
   git commit -m "Add feature: description"
   git push origin feature/your-feature-name
   ```

6. **Create a Pull Request**
   - Describe what your changes do
   - Reference any related issues
   - Include example output if applicable

### Code of Conduct

- Be respectful and inclusive
- Provide constructive feedback
- Focus on the science and methodology
- Acknowledge contributions from others

---

## Citation

If you use this workflow in your research, please cite:

### Primary Citation

**Manuscript:**
```
García-Estrada, D.A., et al. (2026). When Fungi Speak the Language of Plants: Shared Phytohormones with Divergent Meanings...
```

**BibTeX:**
```bibtex
@article{garcia2026,
  title={When Fungi Speak the Language of Plants: Shared Phytohormones with Divergent Meanings},
  author={García-Estrada, David Alberto and ...},
  journal={ASM Microbiology Society},
  year={2026},
  note={Manuscript submitted}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

This work was supported by:

- **Funding**:
- **Databases**: UniProt, NCBI, JGI MycoCosm teams
- **Open-Source Community**: Developers of all bioinformatics tools used

We thank:
- The Bioconductor and CRAN communities for R package development
- All contributors who provided feedback and suggestions

Special acknowledgments:
- Reviewers for their constructive feedback
- The open science community for data sharing

---

## Contact

**David Alberto García Estrada**  
Researcher  
Center for Research in Advanced Materials (CIMAV)  
Chihuahua, México

📧 **Email**: david.garcia@cimav.edu.mx  
🔬 **ORCID**: [0009-0007-1169-5329](https://orcid.org/0009-0007-1169-5329)
💼 **ResearchGate**: [David-Garcia-Estrada](https://www.researchgate.net/profile/David-Garcia-Estrada)

**For technical questions:**
- Open an issue on [GitHub Issues](https://github.com/DavidAlberto/Phytohormones-Fungi/issues)
- Check existing issues before creating new ones

**For collaboration inquiries:**
- Contact via email with subject: "Collaboration - Phytohormone Phylogenomics"

**Repository:**
- 🔗 GitHub: https://github.com/DavidAlberto/Phytohormones-Fungi
- 📦 Docker Hub: https://hub.docker.com/repository/docker/davidalbertoge/hormone-analysis/general

---

## FAQ

**Q: Can I use this pipeline for other organisms besides fungi?**  
A: Yes! The pipeline is designed to be taxonomically agnostic. Just replace the fungal protein database with your organisms of interest and adjust metadata accordingly.

**Q: How long does the pipeline take to run?**  
A: For a typical gene family with 100-200 sequences: BLAST (5-30 min), alignment (30 min - 4 hours), IQ-TREE (1-12 hours depending on model complexity). Total: 2-16 hours.

**Q: Do I need a GPU for this pipeline?**  
A: Not for the phylogenetic analysis. GPU is only recommended for AlphaFold structural modeling, which is optional.

**Q: Can I run this on Windows?**  
A: Yes, using Windows Subsystem for Linux (WSL2) or Docker. We recommend Docker for Windows users.

**Q: What if my gene family has >1000 sequences?**  
A: Consider using approximate methods (FastTree instead of IQ-TREE) or representative sampling. IQ-TREE can handle large datasets but may require substantial time.

**Q: Is there a maximum number of sequences I can analyze?**  
A: Practical limits depend on available RAM and time. We've successfully analyzed datasets with 200+ sequences. For >500 sequences, consider filtering or sampling strategies.

---

**Last Updated**: February 2026  
**Maintainer**: David Alberto García Estrada  
**Status**: Active Development

---

<p align="center">
  ⭐ Star this repository if you find it useful! ⭐
</p>
