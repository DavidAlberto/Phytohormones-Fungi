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
