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

# Add ORF names as column to res_df
res_df$ORF <- rownames(res_df)

# Add the mapping to results
res_df <- merge(res_df, gene_map, by = "ORF", all.x = TRUE)

# Defining significant genes 
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()

# Define our background list of genes to compare to
# Rembember that ORA needs an "interesting gene" set and a background set.
all_genes <- res_df %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()

# GO Enrichment 
ego_bp <- enrichGO(gene = sig_genes,
                   universe = all_genes,
                   OrgDb = org.Sc.sgd.db,
                   ont = "BP",
                   keyType = "ORF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = FALSE)
head(as.data.frame(ego_bp))

# KEGG analysis
kegg_enrich <- enrichKEGG(gene = sig_genes,
                          organism = 'sce',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

head(as.data.frame(kegg_enrich))

# Dot Plots
dotplot(ego_bp, showCategory = 20, title = "GO Biological Process")

dotplot(kegg_enrich, showCategory = 15, title = "KEGG Pathway Enrichment")

# Bar plot
barplot(ego_bp, showCategory = 15, title = "GO Biological Process")

barplot(kegg_enrich, showCategory = 15, title = "KEGG Biological Process")

# Comparing Upregulated and Downregulated Genes
upregulated_orf <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()

downregulated_orf <- res_df %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>%
  pull(ORF) %>%
  na.omit() %>%
  unique()


ego_bp_up <- enrichGO(gene = upregulated_orf,
                      universe = all_genes,
                      OrgDb = org.Sc.sgd.db,
                      keyType = "ORF",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

ego_bp_down <- enrichGO(gene = downregulated_orf,
                        universe = all_genes,
                        OrgDb = org.Sc.sgd.db,
                        keyType = "ORF",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)

dotplot(ego_bp_up, showCategory = 15, title = "GO BP - Upregulated Genes")

dotplot(ego_bp_down, showCategory = 15, title = "GO BP - Downregulated Genes")

kegg_up <- enrichKEGG(gene = upregulated_orf,
                      universe = all_genes,
                      organism = "sce",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

dotplot(kegg_up, showCategory = 15, 
        title = "KEGG Pathways - Upregulated Genes")

kegg_down <- enrichKEGG(gene = downregulated_orf,
                        universe = all_genes,
                        organism = "sce",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)

dotplot(kegg_down, showCategory = 15, 
        title = "KEGG Pathways - Downregulated Genes")

gene_list_go <- list(
  Upregulated = upregulated_orf,
  Downregulated = downregulated_orf
)

compare_go <- compareCluster(geneCluster = gene_list_go,
                             fun = "enrichGO",
                             OrgDb = org.Sc.sgd.db,
                             keyType = "ORF",
                             ont = "BP",
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05)

dotplot(compare_go, showCategory = 10, x = "Cluster") +
  labs(title = "GO BP Enrichment Comparison (Up vs Down)")

gene_list_kegg <- list(
  Upregulated = upregulated_orf,
  Downregulated = downregulated_orf
)

compare_kegg <- compareCluster(geneCluster = gene_list_kegg,
                               fun = "enrichKEGG",
                               organism = "sce",
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.05)

dotplot(compare_kegg, showCategory = 10, x = "Cluster") +
  labs(title = "KEGG Pathway Enrichment Comparison (Up vs Down)")
