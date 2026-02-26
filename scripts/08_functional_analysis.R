# Perform GO, KEGG enrichment (ORA) and GSEA
# using Wald results (Stage_3 vs Stage_1)

library(clusterProfiler)
library(enrichplot)
library(org.Sc.sgd.db)
library(DOSE)
library(tidyverse)
library(here)

# Load in figures directory
figures_dir <- here("figures")

# load in results directory
results_dir <- here("results")

# Load Wald results (Stage_3 vs Stage_1)
# Repeat for each stage comparison
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

png(file.path(figures_dir, "GOenrich_DP.png"),
    width = 1000,
    height = 800,
    res = 150)

dotplot(ego_bp, showCategory = 20,
        title = "GO Biological Process - Stage 3 vs Stage 1")

dev.off()

# KEGG analysis
kegg_enrich <- enrichKEGG(gene = sig_genes,
                          organism = 'sce',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

png(file.path(figures_dir, "KEGGenrich_DP.png"),
    width = 1000,
    height = 800,
    res = 150)

dotplot(kegg_enrich, showCategory = 15,
        title = "KEGG Pathway Enrichment - Stage 3 vs Stage 1")

dev.off()

# Bar plot
png(file.path(figures_dir, "GOenrich_BP.png"),
    width = 1000,
    height = 800,
    res = 150)

barplot(ego_bp, showCategory = 15, title = "GO Biological Process - Stage 3 vs Stage 1")

dev.off()

png(file.path(figures_dir, "KEGGenrich_BP.png"),
    width = 1000,
    height = 800,
    res = 150)

barplot(kegg_enrich, showCategory = 15, title = "KEGG Pathway Enrichment - Stage 3 vs Stage 1")
dev.off()

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

png(file.path(figures_dir, "GOenrich_upvsdown.png"),
    width = 1000,
    height = 800,
    res = 150)

dotplot(compare_go, showCategory = 7) +
  labs(title = "GO Biological Process (Up vs Down)")

dev.off()

# KEGG Enrichment Comparison (Up vs Down)
compare_kegg <- compareCluster(
  geneCluster = gene_sets,
  fun = "enrichKEGG",
  organism = "sce",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

png(file.path(figures_dir, "KEGGenrich_upvsdown.png"),
    width = 1000,
    height = 800,
    res = 150)

dotplot(compare_kegg, showCategory = 7) +
  labs(title = "KEGG Pathway Enrichment (Up vs Down)")

dev.off()

# Repeat the same with the clusters from the LRT

# GSEA (Ranked Gene List)
# Repeat this with LRT CLusters 
res_clean <- res_df

# Remove genes with NA test statistics 
res_clean <- res_clean[!is.na(res_clean$stat), ]

# Remove duplicated ORF entries to ensure one statistice per gene
res_clean <- res_clean[!duplicated(res_clean$ORF), ]

# Create ranked numeric vector of test statistics
gene_list <- res_clean$stat

# gene identifiers which are the ORF names
names(gene_list) <- res_clean$ORF

# clean any NA's
gene_list <- gene_list[!is.na(gene_list)]

# Sort genes in decreasing order so that upregulated genes are at the top
# and downregulated genes are at the bottom
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_go <- gseGO(geneList = gene_list,
                 OrgDb = org.Sc.sgd.db,
                 keyType = "ORF",
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

write.csv(as.data.frame(gsea_kegg),
          file = file.path(results_dir, "GSEA_GO_results.csv"),
          row.names = FALSE)

png(file.path(figures_dir, "GSEA_GODP.png"),
    width = 1000,
    height = 800,
    res = 150)

dotplot(gsea_go, showCategory = 15,
        title = "GSEA - GO Biological Process")

dev.off()

p <- gseaplot2(gsea_go, geneSetID = 1)
ggsave(
  filename = file.path(figures_dir, "GSEA_GO_geneSet1.png"),
  plot = p,
  width = 8,
  height = 6
)

gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = "sce",
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

write.csv(as.data.frame(gsea_kegg),
          file = file.path(results_dir, "GSEA_KEGG_results.csv"),
          row.names = FALSE)

png(file.path(figures_dir, "GSEA_KEGGDP.png"),
    width = 1000,
    height = 800,
    res = 150)

dotplot(gsea_kegg, showCategory = 15,
        title = "GSEA - KEGG Pathways")
dev.off()

p <- gseaplot2(gsea_kegg, geneSetID = 1)

ggsave(
  filename = file.path(figures_dir, "GSEA_KEGG_geneSet1.png"),
  plot = p,
  width = 8,
  height = 6
)
as.data.frame(gsea_kegg) %>% head()
