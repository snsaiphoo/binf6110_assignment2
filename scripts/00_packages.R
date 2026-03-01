# CRAN packages
cran_packages <- c(
  "tidyverse",
  "here",
  "pheatmap"
)

installed_cran <- cran_packages %in% rownames(installed.packages())

if (any(!installed_cran)) {
  install.packages(cran_packages[!installed_cran])
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "tximport",
  "DESeq2",
  "GenomicFeatures",
  "org.Sc.sgd.db",
  "AnnotationDbi",
  "DEGreport",
  "apeglm",
  "clusterProfiler",
  "enrichplot",
  "DOSE"
)

installed_bioc <- bioc_packages %in% rownames(installed.packages())

if (any(!installed_bioc)) {
  BiocManager::install(bioc_packages[!installed_bioc], ask = FALSE)
}
