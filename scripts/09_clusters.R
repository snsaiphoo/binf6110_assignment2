# Perform differential expression analysis (LRT + Wald),
# generate PCA, MA, volcano plots, and heatmaps.

library(DESeq2)
library(org.Sc.sgd.db)
library(AnnotationDbi)
library(tidyverse)
library(here)
library(DEGreport)
library(apeglm)
library(pheatmap)
library(enrichplot)
library(clusterProfiler)


# Load location of results directory
results_dir <- here("results")

# Create folder for figures, and define location
figures_dir <- here("figures")
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir)
}

# Load DESeq2 dataset
dds <- readRDS(here("dds_object.rds"))

# Pre-filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq model 
dds <- DESeq(dds)

# Variance-stabilizing transformation (for visualization)
vsd <- vst(dds, blind = FALSE)

# Extract matrix
vsd_mat <- assay(vsd)

########## Likelihood ratio test (LRT) ######################
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)

res_LRT <- results(dds_lrt)

res_LRT_tb <- as.data.frame(res_LRT) %>%
  rownames_to_column("gene") %>%
  as_tibble()

sigLRT_genes <- res_LRT_tb %>%
  filter(padj < 0.05)

nrow(sigLRT_genes)

# Cluster significant LRT genes 
clustering_sig_genes <- sigLRT_genes %>%
  arrange(padj) %>%
  head(n = 1000)

rld_mat <- assay(vsd)

cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]

meta <- as.data.frame(colData(dds))

clusters <- degPatterns(cluster_rlog,
                        metadata = meta,
                        time = "condition",
                        col = NULL)

cluster_assignments <- clusters$df %>%
  dplyr::select(genes, cluster) %>%
  distinct()

unique_clusters <- unique(cluster_assignments$cluster)

# Heatmaps

for (cl in unique_clusters) {
  
  # Get genes in this cluster
  cluster_genes <- cluster_assignments %>%
    filter(cluster == cl) %>%
    pull(genes)
  
  # Convert dots back to dashes
  cluster_genes <- gsub("\\.", "-", cluster_genes)
  
  # Subset transformed matrix
  # Convert dots back to dashes
  mat <- vsd_mat[cluster_genes, , drop = FALSE]
  
  # Heatmap
  png(file.path(figures_dir,
                paste0("Heatmap_LRT_cluster_", cl, ".png")),
      width = 1000,
      height = 1000,
      res = 150)
  
  pheatmap(mat,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = meta,
           show_rownames = FALSE,
           main = paste("LRT Cluster", cl)
  )
  
  dev.off()
}

## Functional Analysis 

cluster1_genes <- cluster_assignments$genes[cluster_assignments$cluster == 1]
cluster2_genes <- cluster_assignments$genes[cluster_assignments$cluster == 2]
cluster3_genes <- cluster_assignments$genes[cluster_assignments$cluster == 3]
cluster4_genes <- cluster_assignments$genes[cluster_assignments$cluster == 4]

# GO Enrichment
for (i in 1:4) {
  
  go_obj <- get(paste0("cluster", i,"_genes"))
  
  # Fix dot to dash if needed
  cluster_genes <- gsub("\\.", "-", go_obj)
  
  if (length(cluster_genes) < 5) {
    message("Cluster ", i, " has too few genes — skipping")
    next
  }
  
  # Run GO enrichment
  go_res <- enrichGO(
    gene          = cluster_genes,
    OrgDb         = org.Sc.sgd.db,
    keyType       = "ORF",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )

  # Skip if no enrichment
  if (nrow(as.data.frame(go_res)) == 0) {
    next
  }
  
  # Save plot
  png(file.path(figures_dir,
                paste0("GOenrich_LRT_cluster_", i, ".png")),
      width = 1000,
      height = 800,
      res = 150)
  
  print(
  dotplot(go_res) +
      labs(title = paste("GO BP Enrichment - LRT Cluster", i)))
  
  dev.off()
}

## KEGG LRT CLUSTERS #####

kegg_cluster1 <- enrichKEGG(
  gene         = cluster1_genes,
  organism     = "sce",
  pvalueCutoff = 0.05
)

kegg_cluster2 <- enrichKEGG(
  gene         = cluster2_genes,
  organism     = "sce",
  pvalueCutoff = 0.05
)

kegg_cluster3 <- enrichKEGG(
  gene         = cluster3_genes,
  organism     = "sce",
  pvalueCutoff = 0.05
)

kegg_cluster4 <- enrichKEGG(
  gene         = cluster4_genes,
  organism     = "sce",
  pvalueCutoff = 0.05
)

for (i in 1:4) {
  
  kegg_obj <- get(paste0("kegg_cluster", i))
  
  # If no enrichment found, skip
  if (is.null(kegg_obj) || nrow(as.data.frame(kegg_obj)) == 0) {
    next
  }
  
  png(file.path(figures_dir, paste0("KEGGenrich_LRT_cluster", i, ".png")),
      width = 1000,
      height = 800,
      res = 150)
  
  print(dotplot(kegg_obj))
  
  dev.off()
}

########## GO GSEA FOR LRT CLUSTERS ######################

# Use LRT results object
gene_list <- res_LRT$stat
names(gene_list) <- rownames(res_LRT)  

# Remove NA
gene_list <- gene_list[!is.na(gene_list)]

# Sort decreasing
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_go <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Sc.sgd.db,
  keyType       = "ORF",
  ont           = "BP",
  minGSSize     = 10,
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

png(file.path(figures_dir, "GSEA_GO_LRT.png"),
    width = 1200,
    height = 900,
    res = 150)

print(dotplot(gsea_go, showCategory = 15))

dev.off()

png(file.path(figures_dir, "GSEA_GO_LRT_top_pathway.png"),
    width = 1200,
    height = 900,
    res = 150)

print(
  gseaplot2(
    gsea_go,
    geneSetID = 1,
    title = gsea_go@result$Description[1]
  )
)

dev.off()


### KEGG GSEA for LRT-ranked gene list ####

gene_list <- res_LRT$stat
names(gene_list) <- rownames(res_LRT)   # OR res_lrt$ORF
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_kegg <- gseKEGG(
  geneList      = gene_list,
  organism      = "sce",
  minGSSize     = 10,
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

png(file.path(figures_dir, "GSEA_KEGG_LRT_dotplot.png"),
    width = 1200,
    height = 900,
    res = 150)

dev.off()

png(file.path(figures_dir, "GSEA_KEGG_LRT_top_pathway.png"),
    width = 1200,
    height = 900,
    res = 150)

print(
  gseaplot2(
    gsea_kegg,
    geneSetID = 1,
    title = gsea_kegg@result$Description[1]
  ))


dev.off()
