# Perform differential expression analysis (LRT + Wald),
# generate PCA, MA, volcano plots, and heatmaps.

library(DESeq2)
library(tidyverse)
library(here)
library(DEGreport)
library(apeglm)
library(pheatmap)

# Create results directory if it doesn't exist
results_dir <- here("results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Load DESeq2 dataset
dds <- readRDS(here("dds_object.rds"))

# Pre-filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq model 
dds <- DESeq(dds)

# Variance-stabilizing transformation (for visualization)
vsd <- vst(dds, blind = FALSE)

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

clusters <- degPatterns(cluster_rlog,
                        metadata = meta,
                        time = "condition",
                        col = NULL)


############# Wald Tests #############

comparisons <- list(
  c("condition","Stage_2","Stage_1"),
  c("condition","Stage_3","Stage_1"),
  c("condition","Stage_3","Stage_2")
)

# PCA Plots
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA Plot of Samples") +
  theme_minimal()

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
  plotMA(resLFC,
         ylim = c(-5,5),
         main = paste("MA Plot:", contrast_name))
  
  # Volcano Plot
  res_df <- as.data.frame(resLFC)
  
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
    ggplot(res_df,
           aes(log2FoldChange,
               -log10(padj),
               color = significant)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("Down"="blue",
                                    "Not Sig"="gray",
                                    "Up"="red")) +
      labs(title = paste("Volcano Plot:", contrast_name),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value") +
      theme_minimal()
  )
  
  cat("Significant genes in", contrast_name, ":",
      sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
}


