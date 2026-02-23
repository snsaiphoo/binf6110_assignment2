library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(here)

#BiocManager::install(c("tximport","DESeq2"))
#BiocManager::install("GenomicFeatures")
library(GenomicFeatures)

salmon_dir <- here("salmon_output")
salmon_dir

sample_table <- data.frame(
  sample = c("SRR10551657","SRR10551658","SRR10551659",
             "SRR10551660","SRR10551661","SRR10551662",
             "SRR10551663","SRR10551664","SRR10551665"),
  condition = c("Stage_3","Stage_3","Stage_3",
                "Stage_2","Stage_2","Stage_2",
                "Stage_1","Stage_1","Stage_1")
)

# It replaces the row numbers with the sample IDs.
rownames(sample_table) <- sample_table$sample

files <- file.path(salmon_dir,
                   sample_table$sample,
                   "quant.sf")

names(files) <- sample_table$sample

file.exists(files)

# Obtain transcript file

gtf_file <- here("data", "genomic.gtf")
gtf_file

txdb <- makeTxDbFromGFF(gtf_file)

k <- keys(txdb, keytype = "TXNAME")

tx2gene <- select(txdb,
                  keys = k,
                  columns = "GENEID",
                  keytype = "TXNAME")

head(tx2gene)
dim(tx2gene)

# To fix for naming inconsistencies
tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)
head(tx2gene)

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

dds <- dds[rowSums(counts(dds)) > 10, ]

dim(dds)

dds <- DESeq(dds)

resultsNames(dds)

#Generate a results table for each stage
res_2_vs_1 <- results(dds,
                      contrast = c("condition",
                                   "Stage_2",
                                   "Stage_1"))
res_3_vs_1 <- results(dds,
                      contrast = c("condition",
                                   "Stage_3",
                                   "Stage_1"))
res_3_vs_2 <- results(dds,
                      contrast = c("condition",
                                   "Stage_3",
                                   "Stage_2"))

vsd <- vst(dds)
plotPCA(vsd, intgroup = "condition")


subset(res_2_vs_1, padj<.05)
plotMA(res_2_vs_1, ylim = c(-5,5))
plotMA(res_3_vs_1, ylim = c(-5,5))

dim(dds)
resultsNames(dds)
summary(res_3_vs_1)







