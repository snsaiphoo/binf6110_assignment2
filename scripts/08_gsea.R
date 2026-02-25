library(BiocManager)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(org.Sc.sgd.db)
library(ggplot2)
library(DOSE)
library(tidyverse)


res_clean <- as.data.frame(resLFC)

res_clean$ORF <- rownames(res_clean)

# Remove NA values
res_clean <- res_clean[!is.na(res_clean$log2FoldChange), ]

# Remove duplicates (keep first occurrence)
res_clean <- res_clean[!duplicated(res_clean$ORF), ]

gene_list <- res_clean$stat
names(gene_list) <- res_clean$ORF

gene_list <- gene_list[!is.na(gene_list)]
gene_list <- sort(gene_list, decreasing = TRUE)


gsea_go <- gseGO(geneList = gene_list,
                 OrgDb = org.Sc.sgd.db,
                 keyType = "ORF",
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

dotplot(gsea_go, showCategory = 15,
        title = "GSEA - GO Biological Process")

## KEGG

gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = "sce",
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

dotplot(gsea_kegg, showCategory = 15,
        title = "GSEA - KEGG Pathways")
