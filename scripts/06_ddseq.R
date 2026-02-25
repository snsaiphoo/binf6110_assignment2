library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(here)

#BiocManager::install(c("tximport","DESeq2"))
#BiocManager::install("GenomicFeatures")
library(GenomicFeatures)

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

# Obtain transcript file

gtf_file <- here("data", "genomic.gtf")
gtf_file

txdb <- makeTxDbFromGFF(gtf_file)

k <- keys(txdb, keytype = "TXNAME")

tx2gene <- select(txdb,
                  keys = k,
                  columns = "GENEID",
                  keytype = "TXNAME")

head(tx2gene)
dim(tx2gene)

# To fix for naming inconsistencies
tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)
head(tx2gene)

txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

dim(txi$counts)

# Creating the DESeq2 dataset

dds <- DESeqDataSetFromTximport(
  txi,
  colData = sample_table,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

dim(dds)

dds <- DESeq(dds)
vsd <- vst(dds)

comparisons <- list(
  list(type="coef", name="condition_Stage_2_vs_Stage_1"),
  list(type="coef", name="condition_Stage_3_vs_Stage_1"),
  list(type="contrast", name=c("condition","Stage_3","Stage_2"))
)

# PCA Plots
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

library(ggplot2)

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Samples") +
  theme_minimal()

# Visualizations
library(apeglm)
library(pheatmap)
library(ggplot2)

for (comp in comparisons) {
  
  if (comp$type == "coef") {
    
    cat("\nProcessing:", comp$name, "\n")
    
    resLFC <- lfcShrink(dds,
                        coef = comp$name,
                        type = "apeglm")
    
    plot_title <- comp$name
    
  } else {
    
    cat("\nProcessing: Stage_3 vs Stage_2\n")
    
    resLFC <- results(dds,
                      contrast = comp$name)
    
    plot_title <- "Stage_3_vs_Stage_2"
  }
  
  # MA Plot
  plotMA(resLFC,
         ylim = c(-5,5),
         main = paste("MA Plot:", plot_title))
  
  # Volcano Plot
  res_df <- as.data.frame(resLFC)
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
      labs(title = paste("Volcano Plot:", plot_title),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value") +
      theme_minimal()
  )
  
  # Heatmap (Top 20)
  res_clean <- na.omit(resLFC)
  top_genes <- head(order(res_clean$padj), 20)
  gene_names <- rownames(res_clean)[top_genes]
  
  mat <- assay(vsd)[gene_names, ]
  
  annotation_df <- data.frame(Stage = sample_table$condition)
  rownames(annotation_df) <- rownames(sample_table)
  
  pheatmap(mat,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = annotation_df,
           show_rownames = TRUE,
           show_colnames = FALSE,
           main = paste("Top 20 DE Genes:", plot_title))
}

# Differential Counts
cat("Significant genes:",
    sum(resLFC$padj < 0.05, na.rm=TRUE), "\n")






