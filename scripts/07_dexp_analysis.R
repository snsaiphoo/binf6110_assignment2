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

# Save full LRT results
write.csv(res_LRT_tb,
          file = file.path(results_dir, "LRT_full_results.csv"),
          row.names = FALSE)

# Save significant LRT genes
write.csv(sigLRT_genes,
          file = file.path(results_dir, "LRT_significant_genes.csv"),
          row.names = FALSE)

# Cluster significant LRT genes 
clustering_sig_genes <- sigLRT_genes %>%
  arrange(padj) %>%
  head(n = 1000)

rld_mat <- assay(vsd)

cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]

meta <- as.data.frame(colData(dds))

png(file.path(figures_dir, "LRT.png"),
    width = 1000,
    height = 800,
    res = 150)

clusters <- degPatterns(cluster_rlog,
                        metadata = meta,
                        time = "condition",
                        col = NULL)

dev.off()

############# Wald Tests #############

comparisons <- list(
  c("condition","Stage_2","Stage_1"),
  c("condition","Stage_3","Stage_1"),
  c("condition","Stage_3","Stage_2")
)

# PCA Plots
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
# Get percent variance explained by the top two principal components
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Samples") +
  theme_minimal()

ggsave(filename = file.path(figures_dir, "PCA_plot.png"),
       plot = pca,
       width = 8,
       height = 6 
       )

# Visualizations
for (comp in comparisons) {
  
  contrast_name <- paste(comp[2], "vs", comp[3])
  
  # For Stage 3 vs Stage 1 and Stage 2 vs Stage 1
  if (comp[3] == "Stage_1") {
    
    coef_name <- paste0("condition_", comp[2], "_vs_", comp[3])
    
    resLFC <- lfcShrink(dds,
                        coef = coef_name,
                        type = "apeglm")
    
  } else {
    
    # For Stage_3 vs Stage_2 
    resLFC <- results(dds,
                      contrast = comp)
  }
  
  # MA Plot
  
  png(file.path(figures_dir, paste0("MA ", contrast_name, ".png")),
      width = 1000,
      height = 800,
      res = 150)
  
  plotMA(resLFC,
         ylim = c(-5,5),
         main = paste("MA Plot:", contrast_name))
  
  dev.off()
  
  # Volcano Plot
  res_df <- as.data.frame(resLFC)
  res_df$gene <- rownames(res_df)
  
  # Save full Wald results
  write.csv(res_df,
            file = file.path(results_dir,
                             paste0("Wald_full_", gsub(" ", "_", contrast_name), ".csv")),
            row.names = FALSE)
  
  # Save significant genes only
  sig_df <- res_df %>%
    filter(padj < 0.05,
           abs(log2FoldChange) > 1)
  
  write.csv(sig_df,
            file = file.path(results_dir,
                             paste0("Wald_significant_", gsub(" ", "_", contrast_name), ".csv")),
            row.names = FALSE)
  
  res_df$gene <- rownames(res_df)
  
  res_df$significant <- ifelse(res_df$padj < 0.05 &
                                 abs(res_df$log2FoldChange) > 1,
                               ifelse(res_df$log2FoldChange > 0,
                                      "Up","Down"),
                               "Not Sig")
  
  res_df <- na.omit(res_df)
  print(
    vp <- ggplot(res_df,
           aes(log2FoldChange,
               -log10(padj),
               color = significant)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("Down"="blue",
                                    "Not Sig"="gray",
                                    "Up"="red")) +
      labs(title = paste("Volcano Plot:", contrast_name),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value") + theme_minimal()
    )
  
  ggsave(
    filename = file.path(figures_dir, paste0("Volcano_", contrast_name, ".png")),
    plot = vp,
    width = 8,
    height = 6
  )
  
  # Heatmaps
  # Remove NA padj values first
  res_no_na <- resLFC[!is.na(resLFC$padj), ]
  
  # Select top 15 genes by adjusted p-value
  top_genes <- rownames(res_no_na)[order(res_no_na$padj)][1:15]
  
  # Map ORF (locus) to gene symbol
  gene_symbols <- mapIds(org.Sc.sgd.db,
                         keys = top_genes,
                         column = "GENENAME",
                         keytype = "ORF",
                         multiVals = "first")
  
  # Replace NA symbols with locus if missing
  gene_symbols[is.na(gene_symbols)] <- top_genes[is.na(gene_symbols)]
  
  # Create combined label: Locus_Gene 
  heatmap_labels <- paste0(top_genes, "_", gene_symbols)
  
  # Subset transformed matrix
  mat <- vsd_mat[top_genes, ]
  
  # Replace rownames in matrix
  rownames(mat) <- heatmap_labels
  
  # Optional: ensure annotation rownames match colnames
  annotation_df <- as.data.frame(colData(dds))
  
  # Save heatmap
  png(file.path(figures_dir,
                paste0("Heatmap_", gsub(" ", "_", contrast_name), ".png")),
      width = 1000,
      height = 1000,
      res = 150)
  
  pheatmap(mat,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = annotation_df,
           show_rownames = TRUE,
           show_colnames = FALSE,
           main = paste("Top 20 DE Genes:", contrast_name)
  )
  
  dev.off()
  
  cat("Significant genes in", contrast_name, ":",
      sum(res_df$padj < 0.05, na.rm = TRUE), "\n",
      file = file.path(figures_dir, "sig_genes.txt"),
      append = TRUE)
}


