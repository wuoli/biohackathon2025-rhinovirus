# Codon Usage Bias

This folder contains the codon usage analysis used to investigate patterns of **codon usage bias in human rhinovirus (HRV)** and their potential relationship to immune avoidance and translational selection.

## Folder Structure

* **`Rhinovirus_Genome.R`** — Main analysis of codon frequencies across 16 HRV strains. Calculates average codon usage and variation across strains and compares HRV codon usage with human codon usage.
* **`West_Genome.R`** — Codon usage analysis for 12 West Nile Virus (WNV) strains, providing a comparison with human codon usage and the HRV analysis.
* **`CAI_Calculator.R`** — Calculates the Codon Adaptation Index (CAI) for rhinovirus using human codon usage as the reference.
* **`CAI_Graph.R`** — Summarizes and visualizes CAI values across HRV and WNV strains and performs a Wilcoxon test comparing the two groups.

## Analysis Workflow

The analysis examined viral codon usage through two complementary approaches: **codon frequency** and **Codon Adaptation Index (CAI)**.

`Rhinovirus_Genome.R` calculates codon frequencies across 16 HRV strains and compares them with human codon usage. `West_Genome.R` performs the corresponding analysis for WNV, allowing the two viruses to be compared.

The CAI analysis evaluates how closely viral codon usage aligns with human codon preferences. `CAI_Calculator.R` calculates CAI values using human codon usage as the reference, while `CAI_Graph.R` summarizes the resulting values and compares CAI between HRV and WNV.

Together, these analyses were used to investigate whether HRV codon usage patterns may reflect **immune avoidance pressures in addition to translational selection**.
