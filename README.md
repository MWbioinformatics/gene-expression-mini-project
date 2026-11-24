# Comparative Gene Expression Analysis

## Project Overview
This project explores differential gene expression in human airway smooth muscle cells following dexamethasone treatment, using the Bioconductor `airway` RNA-seq dataset.

The analysis was conducted in R using Bioconductor and the DESeq2 package.

## Objectives
- Explore and assess RNA-seq count data.
- Perform exploratory data analysis and sample-level quality checks.
- Normalise gene expression counts.
- Investigate sample relationships using principal component analysis (PCA).
- Identify differentially expressed genes between treated and untreated samples.
- Visualise and interpret differential expression results.

## Tools and Packages
- R and RStudio
- tidyverse
- Bioconductor
- DESeq2
- SummarizedExperiment
- ggplot2
- org.Hs.eg.db

## Analysis Workflow
1. Data loading and inspection
2. Count extraction and quality assessment
3. Exploratory data analysis
4. Sample-level visualisation
5. Count normalisation and PCA
6. Differential gene expression analysis
7. Gene annotation and results export

## Key Findings
- Identified 4,834 genes with an adjusted p-value below 0.1.
- 2,610 genes were upregulated following dexamethasone treatment.
- 2,224 genes were downregulated.
- PCA suggested separation between treated and untreated samples.
- Annotated significant genes using Ensembl gene identifiers and human gene symbols.

## Repository Structure
- `gene_expression_analysis.R` – Main analysis script.
- `results/` – Exported differential expression results.

## Dataset
The analysis uses the `airway` dataset from Bioconductor, containing RNA-seq data from eight human airway smooth muscle cell samples across four cell lines, with paired dexamethasone-treated and untreated conditions.

## Notes
This is an educational bioinformatics portfolio project. The analysis demonstrates an introductory RNA-seq workflow and differential expression analysis using R and DESeq2.