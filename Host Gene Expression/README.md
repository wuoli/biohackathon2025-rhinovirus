# Host Gene Expression

This folder contains the host gene expression analysis used to investigate **transcriptional responses to human rhinovirus (HRV) and West Nile Virus (WNV) infection**.

## Folder Structure

* **`DEGs_WNV_HRV.Rmd`** — Main analysis file for the differential gene expression analysis, pathway analysis, and the generation of all reported figures. The analysis compares host gene expression responses to HRV and WNV infection relative to controls.
* **`DEG_microarray_rhinovirus_vs_control.csv`** — Differentially expressed gene results for HRV versus control, generated from `DEGs_WNV_HRV.Rmd`.
* **`DESeq2_results_WNV_vs_control.csv`** — Differentially expressed gene results for WNV versus control, generated from `DEGs_WNV_HRV.Rmd`.
* **`pathwayanalysis_grouped.R`** — Code file for alternative options of pathway analysis, ultimately not used for the project.
* **`enrichmen_genes.csv`** — Gene list used for downstream enrichment analysis.
* **`enrichment_results.csv`** — Results from the gene enrichment analysis.
* **`GSE40718_raw_counts_GRCh38.p13_NCBI.tsv.gz`** — Raw RNA-seq gene count matrix used as the input for the WNV DESeq2 analysis. The analysis selects 5 control and 5 WNV samples from this matrix.
* **`GSM*.CEL.gz`** — Raw microarray CEL files used for the HRV host gene expression analysis.
* **`GSM*_transcripts.expr.txt.gz`** — Individual transcript-level expression files corresponding to the WNV/control samples. These files are included in the repository but are not directly read in the shown DESeq2 workflow; the analysis instead uses the corresponding sample columns from `GSE40718_raw_counts_GRCh38.p13_NCBI.tsv.gz`.

## Analysis Workflow

The analysis used publicly available host gene expression datasets to compare transcriptional responses following **HRV and WNV infection**.

For HRV, raw microarray CEL files were normalized using **RMA (Robust Multi-array Average)** and analyzed using **limma** to identify differentially expressed genes between HRV-infected and control samples.

For WNV, the raw RNA-seq count matrix `GSE40718_raw_counts_GRCh38.p13_NCBI.tsv.gz` was used for **DESeq2** analysis, comparing 5 WNV samples with 5 control samples.

`DEGs_WNV_HRV.Rmd` contains the main analysis and generates the differential expression results, GO pathway enrichment, and figures. The resulting HRV and WNV differential expression datasets are saved as `DEG_microarray_rhinovirus_vs_control.csv` and `DESeq2_results_WNV_vs_control.csv`, respectively.

The resulting differentially expressed genes were also explored with other pathway analysis methods such as GSEA in the files `pathwayanalysis_grouped.R`, `enrichmen_genes.csv`, and `enrichment_results.csv`. This analysis was ultimately used as a sanity check and not included in the paper.

Together, these analyses were used to characterize differences in **host transcriptional responses to HRV and WNV infection**, specifically differences in antiviral and innate immune response pathways.
