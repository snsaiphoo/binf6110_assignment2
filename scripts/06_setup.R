# Import Salmon transcript quantifications, summarize to gene-
# level counts using tximport, and construct a DESeq2 dataset.

# Load libraries
library(tximport)
library(DESeq2)
library(tidyverse)
library(here)
library(GenomicFeatures)

# Create results directory if it doesn't exist
results_dir <- here("results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Define directory for the quant files 
salmon_dir <- here("salmon_output")

# Sample metadata dataframe
# Three biological conditions (velum developmental stages)
# with three replicates per stage
sample_table <- data.frame(
  sample = c("SRR10551657","SRR10551658","SRR10551659",
             "SRR10551660","SRR10551661","SRR10551662",
             "SRR10551663","SRR10551664","SRR10551665"),
  condition = c("Stage_3","Stage_3","Stage_3",
                "Stage_2","Stage_2","Stage_2",
                "Stage_1","Stage_1","Stage_1")
)

# Convert condition to factor (levels determined alphabetically)
sample_table$condition <- factor(sample_table$condition)

rownames(sample_table) <- sample_table$sample

# the locations of all the files
files <- file.path(salmon_dir,
                   sample_table$sample,
                   "quant.sf")

names(files) <- sample_table$sample

# verify the files are there
file.exists(files)

# Obtain transcript file from data folder
gtf_file <- here("data", "genomic.gtf")
gtf_file

# Create TxDb object from annotation
txdb <- makeTxDbFromGFF(gtf_file)

# Extract transcript IDs
k <- keys(txdb, keytype = "TXNAME")

# Create transcript-to-gene mapping table
tx2gene <- AnnotationDbi::select(txdb,
                  keys = k,
                  columns = "GENEID",
                  keytype = "TXNAME")

head(tx2gene)
dim(tx2gene)

# To fix for naming inconsistencies
tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)
head(tx2gene)

# tximport summarizes transcript-level estimates to gene-level counts
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

dim(txi$counts)

# Creating the DESeq2 dataset
dds <- DESeqDataSetFromTximport(
  txi,
  colData = sample_table,
  design = ~ condition
)

# Save DESeq2 object for downstream analysis
saveRDS(dds, file = file.path(results_dir, "dds_object.rds"))
