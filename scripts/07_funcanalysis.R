# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("org.Sc.sgd.db")

library(BiocManager)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(org.Sc.sgd.db)
library(ggplot2)
library(DOSE)
library(tidyverse)

# Convert to dataframe
res_df <- as.data.frame(resLFC)
res_df

orf_ids <- rownames(res_df)

# Map to Entrez IDs
gene_map <- bitr(orf_ids, 
                 fromType = "ORF", 
                 toType = c("ENTREZID", "COMMON"),
                 OrgDb = org.Sc.sgd.db)


