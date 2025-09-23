# *Phytohormones in Filamentous Fungi: Biosynthesis, Perception, and Growth Regulation Response - Pipeline*

This repository contains a reproducible, documented bioinformatics pipeline intended to support the review
“*Phytohormones in Filamentous Fungi: Biosynthesis, Perception, and Growth Regulation Response*”.
The pipeline is designed to help researchers:

- find homologs of phytohormone-related genes across fungal genomes/proteomes,
- retrieve and curate sequences,
- build high-quality multiple sequence alignments (MSAs),
- trim alignments for phylogenetic inference,
- reconstruct robust gene trees,
- generate high-quality structural models with AlphaFold (latest available release),
- and visualize models with ChimeraX.

Recommended tools (examples used in this README): **Blast** for homology searches, **MAFFT** for **MSA**,
**trimAl** for automated trimming, **IQ-TREE** for phylogenetic inference, and **AlphaFold** for protein structure prediction.
> Notes:
> - "This repository is a *complementary resource to the review article: it stores code,*
> *example inputs/outputs, and a documented workflow so readers can reproduce and adapt the analyses.*"
> - "Structural visualization is performed in **ChimeraX**, which should be installed separately".

## Table of Contents

- [Installation](#installation)
  - [Conda environment](#conda-environment)
  - [Docker image](#docker-image)
- [Data requirements](#data-requirements)
- [Pipeline steps](#pipeline-steps)
  - [1. Homology search with Blast](#1-homology-search-with-blast)
  - [2. Sequence retrieval and curation](#2-sequence-retrieval-and-curation)
  - [3. Multiple sequence alignment (MSA)](#3-multiple-sequence-alignment-msa)
  - [4. Alignment trimming](#4-alignment-trimming)
  - [5. Phylogenetic inference](#5-phylogenetic-inference)
  - [6. Protein structure prediction with AlphaFold](#6-protein-structure-prediction-with-alphafold)
  - [7. Structural visualization with ChimeraX](#7-structural-visualization-with-chimerax)
- [Example usage](#example-usage)
- [Requirments](#requirements)
- [Contributing](#contributing)
- [License](#license) 
- [How to cite](#how-to-cite)
- [Acknowledgements](#acknowledgements)
- [Contact](#contact)

## Installation

### Conda environment

Create an enviroment fromt the provided `environment.yml` file:

```bash
conda create --name hormone-analysis --file environment.yml
conda activate hormone-analysis
```

### Docker image

We provide a sample Docker image with all dependencies installed.

To pull the image from DockerHub:

```bash
docker pull davidalbertoge/hormone-analysis:latest
docker run -v $(pwd):/home/ -it davidalbertoge/hormone-analysis:latest
```

## Data requirements

- **Sequence databases**: homology searches we use [UniProt](https://www.uniprot.org), [NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and [JGI MycoCosm fungal genomes](https://mycocosm.jgi.doe.gov/mycocosm/home) databases.
- **AlphaFold**: running in website mode requires internet access to download model parameters and databases ([AlphaFold](https://alphafoldserver.com/)).
- **ChimeraX**: visualization is a separate install([ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html)).

## Pipeline steps

Below are minimal working examples for each major step. Adapt parameters for your dataset.

> Convention: `$PWD` is repository root, `data/raw/` contains input files, `data/results/` will contain outputs.

### 1. Homology search with Blast

Use `blastp` for protein queries against a protein database. Blast analysis you can be performed against online databases or locally.

Example command for local search:

```bash
makeblastdb -in data/raw/fungal_proteomes.fasta -dbtype prot -out fungal_db
blastp -query data/raw/query_protein.fasta -db fungal_db -outfmt 6 -evalue 1e-5 -out data/results/blast_results.tsv
```

For online search, use NCBI's [BLAST web service](https://blast.ncbi.nlm.nih.gov/Blast.cgi), UniProt's [BLAST tool](https://www.uniprot.org/blast/) or JGI MycoCosm's [BLAST tool](https://mycocosm.jgi.doe.gov/mycocosm/blast.jsf).

### 2. Sequence retrieval and curation

Two common tools workflows:

**Local multifasta selection:**

```bash
seqkit grep -f data/raw/sequence_ids.txt data/raw/fungal_proteomes.fasta > data/results/selected_sequences.fasta
```

**NCBI Entrez retrieval:**

```bash
esearch -db protein -query "ACCESSION_OR_QUERY" | efetch -format fasta > data/results/retrieved_sequences.fasta
```

### 3. Multiple sequence alignment (MSA)

Use **MAFFT** for high-quality alignments:

```bash
mafft --auto data/results/selected_sequences.fasta > data/results/alignment.fasta
```

MAFFT provides accurate and scalable methods for protein alignments.

### 4. Alignment trimming

Remove poorly aligned regions with **trimAl**:

```bash
trimal -in data/results/alignment.fasta -out data/results/alignment_trimmed.fasta -automated1
```

`automated1` is a good starting point for trimming.

### 5. Phylogenetic inference

Model selection, ML tree and branch support with **IQ-TREE**:

```bash
iqtree -s data/results/alignment_trimmed.fasta -m MFP -bb 1000 -alrt 1000 -nt AUTO
```

This runs model selection (MFP), 1000 ultrafast bootstrap replicates, and 1000 SH-aLRT tests, for robust support values.

### 6. Protein structure prediction with AlphaFold

**AlphaFold** run in website mode (requires internet access), and outputs PDB files. 

### 7. Structural visualization with ChimeraX

Open ChimeraX and load the PDB file.

## Example usage

Example scripts and input files are provided in the `scripts/` and `data/raw/` directories.

Run each step sequentially, adapting paths and parameters as needed.

## Requirements

- **Compute**: we recommend a machine with at least 16GB RAM and multiple CPU cores for efficient processing.
- **Storage**: ensure sufficient disk space for sequence databases and intermediate files.
- **Internet**: required for online database searches and AlphaFold in website mode.
- **Visualization**: ChimeraX must be installed separately for structural visualization.

## Contributing

Contributions are welcome! Scrips, bugs reports, feature requests, and improvements can be submitted via pull requests or issues. 
Please follow standard GitHub etiquette: clear descriptions, code comments, and testing.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 

## How to cite

If you use this pipeline in your research, please cite the associated review article:
> Author(s). (Year). *Phytohormones in Filamentous Fungi: Biosynthesis, Perception, and Growth Regulation Response*. Journal Name, Volume(Issue), Page numbers. DOI

## Acknowledgements

We thank the open-source community for the development of the tools used in this pipeline.

## Contact

For questions, bugs reports, or suggestions, please open an issue or contact:

- Maintainer: [david.garcia@cimav.edu.mx](david.garcia@cimav.edu.mx)
- Github issues: [GitHub Issues](https://github.com/DavidAlberto/Phytohormones-Fungi/issues)
- Github DavidAlberto: [GitHub Profile](https://github.com/DavidAlberto)

Thank you for using this pipeline! We hope it aids your research into phytohormones in filamentous fungi.
