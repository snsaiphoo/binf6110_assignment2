# #########################
# Version Numbers:
# R 4.5.1
# tidyr 1.3.2  
# here 1.0.2 
# pheatmap 1.0.13 
# ggplot2 4.0.2 
# dplyr 1.2.0
# tximport 1.36.1  
# tidyverse 2.0.0   
# DESeq2 1.48.2
# GenomicFeatures 1.60.0
# org.Sc.sgd.db 3.21.0
# AnnotationDbi 1.70.0 
# DEGreport 1.44.0 
# apeglm 1.30.0  
# clusterProfiler 4.16.0
# enrichplot 1.28.4
# DOSE 4.2.0 
##########################

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
