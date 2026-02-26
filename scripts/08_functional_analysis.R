# Perform GO, KEGG enrichment (ORA) and GSEA
# using Wald results (Stage_3 vs Stage_1)

library(clusterProfiler)
library(enrichplot)
library(org.Sc.sgd.db)
library(DOSE)
library(tidyverse)
library(here)

#Load Wald results (Stage_3 vs Stage_1)
res_df <- read.csv(
  here("results", "Wald_full_Stage_3_vs_Stage_1.csv")
)

# Ensure ORF column exists
res_df$ORF <- res_df$gene

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


dotplot(ego_bp, showCategory = 20,
        title = "GO Biological Process - Stage 3 vs Stage 1")

# KEGG analysis
kegg_enrich <- enrichKEGG(gene = sig_genes,
                          organism = 'sce',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

dotplot(kegg_enrich, showCategory = 15,
        title = "KEGG Pathway Enrichment - Stage 3 vs Stage 1")

# Bar plot
barplot(ego_bp, showCategory = 15, title = "GO Biological Process - Stage 3 vs Stage 1")

barplot(kegg_enrich, showCategory = 15, title = "KEGG Pathway Enrichment - Stage 3 vs Stage 1")

# Comparing Upregulated and Downregulated Genes
gene_sets <- list(
  Upregulated = res_df %>%
    filter(padj < 0.05 & log2FoldChange > 1) %>%
    pull(ORF) %>%
    na.omit() %>%
    unique(),
  
  Downregulated = res_df %>%
    filter(padj < 0.05 & log2FoldChange < -1) %>%
    pull(ORF) %>%
    na.omit() %>%
    unique()
)

# GO Enrichment Comparison (Up vs Down)
compare_go <- compareCluster(
  geneCluster = gene_sets,
  fun = "enrichGO",
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

dotplot(compare_go, showCategory = 7) +
  labs(title = "GO Biological Process (Up vs Down)")

# KEGG Enrichment Comparison (Up vs Down)
compare_kegg <- compareCluster(
  geneCluster = gene_sets,
  fun = "enrichKEGG",
  organism = "sce",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

dotplot(compare_kegg, showCategory = 7) +
  labs(title = "KEGG Pathway Enrichment (Up vs Down)")


# GSEA (Ranked Gene List)

res_clean <- res_df

res_clean <- res_clean[!is.na(res_clean$stat), ]
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
gseaplot2(gsea_go, geneSetID = 1)

gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = "sce",
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

dotplot(gsea_kegg, showCategory = 15,
        title = "GSEA - KEGG Pathways")
gseaplot2(gsea_kegg, geneSetID = 1)
as.data.frame(gsea_kegg) %>% head()
