# Sequence Divergence

This folder contains the sequence divergence analysis used to investigate **genetic variation, mutation rates, phylogenetic relationships, and sequence-level differences in Rhinovirus and West Nile Virus (WNV)**.

## Folder Structure

* **`convert_hittable.sh`** — Extracts unique WNV sequence accessions from the BLAST hit table and creates a cleaned accession list for downstream metadata retrieval.
* **`get_metadata.sh`** — Retrieves WNV sequence metadata, including geographic location and collection date, from GenBank using NCBI EDirect.
* **`rhinovirus_alignment.R`** — Reads Rhinovirus nucleotide sequences, removes duplicate sequences, performs multiple sequence alignment using DECIPHER, and saves the resulting alignment.
* **`rhinovirus_mutation_analysis.R`** — Performs Rhinovirus sequence alignment, phylogenetic analysis, mutation rate estimation, regional mutation rate analysis, UMAP analysis, and sequence divergence comparison with WNV.
* **`rhinovirus_mutation_rate.R`** — Estimates the Rhinovirus evolutionary rate using BEAST and generates the configuration files and outputs required for Bayesian molecular clock analysis.
* **`rhinovirus_phylogenetic_tree.R`** — Constructs Neighbor-Joining and Maximum Likelihood phylogenetic trees for Rhinovirus sequences and estimates sequence divergence and mutation rates using root-to-tip regression.
* **`westnile_mutation_analysis.R`** — Aligns WNV sequences and estimates the WNV mutation rate from sequence differences relative to an early reference sequence.
* **`Results/`** — Contains the generated phylogenetic trees, UMAP visualization, mutation rate plots, sequence alignments, and intermediate analysis files.

## Results

The `Results` folder contains the outputs generated during the sequence divergence analyses.

* **`Divergence Rate Barplot.png`** — Comparison of estimated sequence divergence rates between Rhinovirus and WNV.
* **`PhylogeneticTree.png`** — Phylogenetic tree generated from the analyzed Rhinovirus sequences.
* **`UMAP.png`** — UMAP visualization of Rhinovirus sequences based on genomic sequence distances.
* **`beast_config.xml`** — BEAST configuration file used for Bayesian molecular clock analysis.
* **`fit_opt.rds`** — Saved optimized Maximum Likelihood phylogenetic model.
* **`mutation_lm_fit.rds`** — Saved linear model used for estimating the Rhinovirus mutation rate.
* **`mutations.txt`** — Mutation counts and associated metadata used for mutation rate estimation.
* **`rhinovirus.nex`** — Rhinovirus sequence alignment converted to NEXUS format for phylogenetic analysis.
* **`tree_NJ.rds`** — Saved Neighbor-Joining phylogenetic tree.
* **`umap.txt`** — UMAP coordinates and associated Rhinovirus sequence metadata.

## Analysis Workflow

The analysis used publicly available viral sequence data to investigate **sequence divergence and evolutionary dynamics of Rhinovirus and WNV**.

For **Rhinovirus**, nucleotide sequences were aligned using DECIPHER and analyzed using Neighbor-Joining and Maximum Likelihood phylogenetic methods. Sequence differences were calculated relative to an early reference sequence, and linear regression was used to estimate the mutation rate in substitutions per site per year. Mutation rates were also calculated across individual genomic regions to identify regions with comparatively higher or lower rates of sequence divergence.

Rhinovirus sequence distances were additionally visualized using **UMAP**, allowing the genomic relationships between viral strains to be examined based on sequence similarity. A Bayesian molecular clock analysis using **BEAST** was also performed to estimate the evolutionary rate from temporally sampled sequences.

For **WNV**, sequence accessions were extracted from BLAST results and metadata were retrieved from GenBank. WNV sequences were aligned using DECIPHER, and sequence differences relative to an early reference sequence were used to estimate the mutation rate.

Finally, estimated sequence divergence rates for **Rhinovirus and WNV** were compared to characterize differences in their evolutionary dynamics.

Together, these analyses were used to characterize **viral sequence diversity, phylogenetic relationships, mutation rates, and genomic regions associated with sequence divergence in Rhinovirus and WNV**.
