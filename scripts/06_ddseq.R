library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(here)

#BiocManager::install(c("tximport","DESeq2"))

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
