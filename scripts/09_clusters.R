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

for (cl in unique_clusters) {
  
  message("Processing Cluster: ", cl)
  
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

res_wald <- read.csv(
  here("results", "Wald_full_Stage_3_vs_Stage_1.csv"),
)

res_wald$ORF <- res_wald$gene

for (cl in unique_clusters) {
  
  message("Volcano for Cluster: ", cl)
  
  cluster_genes <- cluster_assignments %>%
    filter(cluster == cl) %>%
    pull(genes)
  
  # Fix dot → dash
  cluster_genes <- gsub("\\.", "-", cluster_genes)
  
  # Keep only genes that exist in Wald results
  cluster_genes <- intersect(cluster_genes, res_wald$ORF)
  
  if (length(cluster_genes) == 0) next
  
  # Create volcano dataframe
  volcano_df <- res_wald %>%
    filter(!is.na(padj))
  
  volcano_df$cluster_member <- ifelse(
    volcano_df$ORF %in% cluster_genes,
    paste0("Cluster_", cl),
    "Other"
  )
  
  p <- ggplot(volcano_df,
              aes(log2FoldChange,
                  -log10(padj),
                  color = cluster_member)) +
    geom_point(alpha = 0.6) +
    theme_minimal() +
    labs(title = paste("Volcano Plot - LRT Cluster", cl),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value")
  
  ggsave(
    filename = file.path(figures_dir,
                         paste0("Volcano_LRT_cluster_", cl, ".png")),
    plot = p,
    width = 8,
    height = 6
  )
}

## Functional Analysis 

library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)

for (cl in unique_clusters) {
  
  message("GO enrichment for Cluster: ", cl)
  
  # Get cluster genes
  cluster_genes <- cluster_assignments %>%
    filter(cluster == cl) %>%
    pull(genes)
  
  # Fix dot to dash
  cluster_genes <- gsub("\\.", "-", cluster_genes)
  
  # Keep only genes present in Wald results
  cluster_genes <- intersect(cluster_genes, res_wald$ORF)
  
  if (length(cluster_genes) < 5) {
    message("Cluster ", cl, " has too few genes — skipping")
    next
  }
  
  # Subset Wald results to cluster genes
  cluster_res <- res_wald %>%
    filter(ORF %in% cluster_genes)
  
  # Define Up and Down genes
  gene_sets <- list(
    Upregulated = cluster_res %>%
      filter(padj < 0.05 & log2FoldChange > 1) %>%
      pull(ORF) %>%
      na.omit() %>%
      unique(),
    
    Downregulated = cluster_res %>%
      filter(padj < 0.05 & log2FoldChange < -1) %>%
      pull(ORF) %>%
      na.omit() %>%
      unique()
  )
  
  # Skip if no significant genes
  if (length(gene_sets$Upregulated) == 0 &&
      length(gene_sets$Downregulated) == 0) {
    message("No significant genes in Cluster ", cl)
    next
  }
  
  # GO Enrichment
  compare_go <- compareCluster(
    geneCluster = gene_sets,
    fun = "enrichGO",
    OrgDb = org.Sc.sgd.db,
    keyType = "ORF",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  # If no enrichment found, skip
  if (is.null(compare_go) || nrow(as.data.frame(compare_go)) == 0) {
    message("No GO enrichment for Cluster ", cl)
    next
  }
  
  # Save plot
  png(file.path(figures_dir,
                paste0("GOenrich_LRT_cluster_", cl, ".png")),
      width = 1000,
      height = 800,
      res = 150)
  
  print(
    dotplot(compare_go, showCategory = 10) +
      labs(title = paste("GO BP Enrichment - LRT Cluster", cl))
  )
  
  dev.off()
  
  # Save table
  write.csv(as.data.frame(compare_go),
            file = file.path(results_dir,
                             paste0("GO_LRT_cluster_", cl, ".csv")),
            row.names = FALSE)
}

########## KEGG GSEA FOR LRT CLUSTERS ######################

library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(dplyr)
library(here)

# Load Wald results (choose which contrast to rank by)
res_wald <- read.csv(
  here("results", "Wald_full_Stage_3_vs_Stage_1.csv")
)

# Keep gene + stat only
res_wald <- res_wald %>%
  dplyr::select(gene, log2FoldChange)

# Remove NA stats
res_wald <- res_wald[!is.na(res_wald$log2FoldChange), ]

# Map ORF -> ENTREZ (KEGG requires ENTREZ IDs)
gene_map <- bitr(res_wald$gene,
                 fromType = "ORF",
                 toType = "ENTREZID",
                 OrgDb = org.Sc.sgd.db)

# Merge mapping
res_wald <- merge(res_wald,
                  gene_map,
                  by.x = "gene",
                  by.y = "ORF")

# Get unique clusters
unique_clusters <- unique(cluster_assignments$cluster)

for (cl in unique_clusters) {
  
  message("Running KEGG GSEA for Cluster: ", cl)
  
  # Get genes in this cluster
  cluster_genes <- cluster_assignments %>%
    filter(cluster == cl) %>%
    pull(genes)
  
  # Subset to cluster genes
  cluster_res <- res_wald %>%
    filter(gene %in% cluster_genes)
  
  # Skip if too few genes
  if (nrow(cluster_res) < 10) {
    message("Cluster ", cl, " has too few genes for GSEA")
    next
  }
  
  # Create ranked vector (ENTREZ IDs required)
  gene_list <- cluster_res$log2FoldChange
  names(gene_list) <- cluster_res$ENTREZID
  
  # Remove missing ENTREZ
  gene_list <- gene_list[!is.na(names(gene_list))]
  
  # Sort decreasing (important!)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Run KEGG GSEA
  gsea_kegg <- gseKEGG(
    geneList = gene_list,
    organism = "sce",
    pvalueCutoff = 0.1,
    minGSSize = 5,
    verbose = FALSE
  )
  
  # If no enrichment found
  if (is.null(gsea_kegg) ||
      nrow(as.data.frame(gsea_kegg)) == 0) {
    
    message("No KEGG GSEA enrichment for Cluster ", cl)
    next
  }
  
  # Save results
  write.csv(as.data.frame(gsea_kegg),
            file = file.path(results_dir,
                             paste0("LRT_GSEA_KEGG_cluster_", cl, ".csv")),
            row.names = FALSE)
  
  # Dotplot
  png(file.path(figures_dir,
                paste0("LRT_GSEA_KEGG_cluster_", cl, ".png")),
      width = 1000,
      height = 800,
      res = 150)
  
  print(dotplot(gsea_kegg,
                showCategory = 10,
                title = paste("KEGG GSEA - LRT Cluster", cl)))
  
  dev.off()
  
  # Enrichment plot of first pathway
  p <- gseaplot2(gsea_kegg, geneSetID = 1)
  
  ggsave(
    filename = file.path(figures_dir,
                         paste0("LRT_GSEA_KEGG_cluster_", cl, "_geneSet1.png")),
    plot = p,
    width = 8,
    height = 6
  )
}
