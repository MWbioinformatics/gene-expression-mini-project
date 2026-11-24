
# ============================================================
# Comparative Gene Expression Analysis
# Dataset: airway (Bioconductor)
# Author: Marc
# Description: Exploratory analysis of dexamethasone-treated
# and untreated human airway smooth muscle cells.
# ============================================================


# DAY 1: SETUP AND DATA LOADING
# ============================================================

# Load required packages
library(tidyverse)
library(airway)
library(DESeq2)
library(SummarizedExperiment)

# Load dataset
data("airway")

# Inspect dataset
airway

# Check dimensions
nrow(airway)  # Number of genes
ncol(airway)  # Number of samples

# Inspect gene and sample identifiers
head(rownames(airway))
head(colnames(airway))

# Inspect sample metadata
colData(airway)


# DAY 2: EXTRACT AND INSPECT COUNTS
# ============================================================

# Extract raw gene expression counts
counts <- assay(airway)

# Inspect count matrix
head(counts)

# Check dimensions and data types
dim(counts)
class(counts)
typeof(counts)

# Check data quality
sum(is.na(counts))
any(is.na(counts))
sum(duplicated(rownames(counts)))

# Summarise count distribution
summary(as.vector(counts))


# DAY 3: EXPLORATORY DATA ANALYSIS
# ============================================================

# Examine treatment groups
table(colData(airway)$dex)

# Examine cell lines
table(colData(airway)$cell)

# Inspect first four samples
head(counts[, 1:4])

# Calculate total counts per sample
# Provides an indication of assigned sequencing reads
colSums(counts)

# Visualise total counts per sample
barplot(
  colSums(counts) / 1e6,
  las = 2,
  cex.names = 0.7,
  ylab = "Total counts (millions)",
  main = "Sequencing depth per sample"
)

# Compare count distributions across samples
# Log2 transformation makes highly skewed counts easier to visualise
# Adding 1 prevents taking the logarithm of zero
boxplot(
  log2(counts + 1),
  las = 2,
  cex.axis = 0.7,
  ylab = "Log2(count + 1)",
  main = "Gene count distributions across samples"
)

# Count detected genes per sample
# TRUE means a gene has at least one assigned read
# colSums counts the TRUE values for each sample
colSums(counts > 0)

# Visualise detected genes per sample
barplot(
  colSums(counts > 0),
  las = 2,
  cex.names = 0.7,
  ylab = "Number of detected genes",
  main = "Genes detected per sample"
)
