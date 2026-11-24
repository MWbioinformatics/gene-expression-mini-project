
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


# DAY 4: DATA VISUALISATION WITH GGPLOT2
# ============================================================

# Create a summary table of sample-level information
sample_qc <- data.frame(
  sample = colnames(counts),
  total_counts = colSums(counts),
  detected_genes = colSums(counts > 0),
  treatment = colData(airway)$dex,
  cell_line = colData(airway)$cell
)

# Inspect the summary table
sample_qc

# Visualise total counts per sample using ggplot2
ggplot(sample_qc, aes(x = sample, y = total_counts / 1e6, fill = treatment)) +
  geom_col() +
  labs(
    title = "Total Assigned Counts per Sample",
    x = "Sample",
    y = "Total counts (millions)",
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Visualise detected genes per sample
ggplot(sample_qc, aes(x = sample, y = detected_genes, fill = treatment)) +
  geom_col() +
  labs(
    title = "Detected Genes per Sample",
    x = "Sample",
    y = "Number of detected genes",
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# DAY 5: NORMALISATION AND PCA
# ============================================================

# Create DESeq2 dataset
dds <- DESeqDataSet(
  airway,
  design = ~ cell + dex
)

# Inspect the DESeq2 object
dds

# Remove genes with very low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Check remaining number of genes
nrow(dds)

# Estimate size factors for library-size normalisation
dds <- estimateSizeFactors(dds)

# Extract normalised counts
normalised_counts <- counts(dds, normalized = TRUE)

# Inspect the normalised count matrix
head(normalised_counts)

# Compare library size factors
sizeFactors(dds)

# Perform variance-stabilising transformation
vsd <- vst(dds, blind = FALSE)

# Generate PCA plot
plotPCA(vsd, intgroup = c("dex", "cell"))
