# Bulk Transcriptomic Analysis of Flor Yeast Velum Development 

## Introduction 
_Saccharomyces cerevisiae_, commonly known as brewer’s yeast, plays a central role in the production of sherry-style wines [1]. During biological wine aging, specialized flor yeast strains form a floating biofilm layer, referred to as a velum, on the surface of fortified wine [1]. This allows the yeast to survive oxidative stress and high ethanol levels that occur during aging. Unlike fermentative yeast strains, which rely primarily on sugars for energy, flor yeast utilizes ethanol and glycerol as carbon sources under oxidative stress, playing a part in both survival and the biochemical characteristics of sherry wines [1].

The dataset analyzed in this study is derived from a transcriptomic investigation of gene expression changes during a 71-day biological wine-aging process [1]. RNA-seq samples were collected at three distinct stages of velum development: Stage 1 (10 days, thin biofilm), Stage 2 (45 days, thin biofilm), and Stage 3 (71 days, mature biofilm). These stages reflect progressive biofilm maturation and physiological adaptation. Assessing gene expression changes across these time points helps identify regulatory mechanisms underlying stress responses, metabolic transitions, and biofilm structural development.

Transcriptomic analysis provides a comprehensive overview of global gene activity and enables investigation of changes in gene expression patterns in response to environmental and physiological stressors [2]. In the context of flor yeast, discovering these transcriptional changes is important for identifying pathways that facilitate oxidative adaptation and prolonged survival during wine aging.

The objective of this analysis is to evaluate differential gene expression across stages and identify biological processes and pathways associated with velum maturation. By using differential expression analysis with functional annotation, this analysis seeks to characterize the transcriptional mechanisms underlying flor yeast adaptation.

### Design Rationale 
Quality control was performed using `FastQC` on raw RNA-seq reads [3] in accordance with established RNA-seq best-practice guidelines [4]. This was to assess per-base sequence quality, GC content distribution, adapter contamination, overrepresented sequences, and duplication levels prior to downstream analysis.

Because transcript isoforms often overlap, short RNA-seq reads may originate from multiple transcripts, making accurate transcript-level quantification essential for reliable downstream analysis [5]. Alignment-free quantification tools such as `Salmon` [6], `Sailfish` [7], and `Kallisto` were evaluated. These methods are computationally efficient and achieve accuracy comparable to or better than that of traditional alignment-based approaches [5].

`Salmon` was selected over `Kallisto` and `Sailfish` because it incorporates more comprehensive bias-correction models (including sequence-specific, GC-content, and positional biases) and a probabilistic two-phase inference framework to improve transcript abundance estimation while maintaining comparable computational efficiency [7]. Additionally, `Salmon` integrates seamlessly with Bioconductor workflows through `tximport`, enabling accurate aggregation of transcript-level estimates to gene-level counts for differential analysis [8, 9].

Differential expression analysis was performed using `DESeq2` [10]. This was selected over `edgeR` [11] and `limma-voom` [12] because it models count data using a negative binomial distribution, which is appropriate for RNA-seq data [10, 13], and performs internal normalization of raw counts, simplifying the workflow. Wald tests were used for pairwise stage comparisons, and likelihood ratio tests (LRTs) were used to detect genes with significant expression changes across multiple stages [10].

To evaluate biological relevance, both over-representation analysis (ORA) and Gene Set Enrichment Analysis (GSEA) were conducted using `clusterProfiler` [14]. ORA identifies enriched Gene Ontology (GO) terms and KEGG pathways based on predefined lists of differentially expressed genes, whereas GSEA assesses enrichment across ranked gene lists without requiring an arbitrary cutoff, allowing detection of coordinated but subtle pathway-level changes [15]. While each statistical and enrichment method has limitations, such as including sensitivity to model assumptions and threshold selection, the combined use of complementary approaches strengthens confidence in the biological interpretation of flor yeast transcriptional regulation during biofilm maturation.

## Methods 
### 1.0 - Data and Tools
#### 1.1 - Containers
All command-line software used in this analysis was executed using containerized environments for reproducibility and version control. Singularity was used to pull pre-built Docker images and convert them into `.sif` container files stored in a dedicated `containers/` directory. 

#### 1.2 - Data Acquisition
The Flor yeast samples were obtained from the NCBI Sequence Read Archive (SRA) and correspond to the yeast biofilm (velum) development study described in [1]. The _Saccharomyces cerevisiae_ samples were collected at three stages of biofilm formation during wine aging:

* **Stage 1** (10 days) – SRR10551665, SRR10551664, SRR10551663
* **Stage 2** (45 days) – SRR10551662, SRR10551661, SRR10551660
* **Stage 3** (71 days) – SRR10551659, SRR10551658, SRR10551657

SRA files were retrieved using the prefetch command from the SRA Toolkit container. The script can be seen in [`01_data.sh`](scripts/01_data.sh).
These SRA files were then converted to FASTQ format using fasterq-dump. The script can be seen in [`02_data.sh`](scripts/02_data.sh). The resulting FASTQ files were compressed for use in downstream quality control and transcript quantification.

#### 1.3 - Reference Genome Reasoning
The RNA-seq data analyzed in this study originate from the _Saccharomyces cerevisiae L-329_ strain described in [1]. However, transcript quantification was performed using the _S288C_ reference transcriptome.

To validate this choice, [`00_refanalysis.sh`](scripts/00_refanalysis.sh) was used to assess sequence conservation between the _L-329_ strain and the _S288C_ reference using BLAST. Results indicated that over 99% of queried sequences aligned to the _S288C_ reference, with an average percent identity of approximately 80%, supporting substantial conservation across coding regions.

Although _L-329_ is the experimental strain, _S288C_ was selected due to its well-annotated and curated genome assembly. Given the high alignment rate and strong sequence conservation, _S288C_ provides an appropriate and biologically relevant reference for transcript quantification and downstream functional enrichment analyses.

The _S288C_ transcript FASTA file was therefore used for downstream analysis.

### 2.0 - Quality Control with FastQC


### 3.0 - Quantification with salmon

### 4.0 - Importing quant files with tximport

### 5.0 - Differential Expression Analysis

#### 5.1 - DESeq2: Wald Test

#### 5.2 - DESeq2: Likelihood Ratio Test (LRT)

### 6.0 - Functional Enrichment Analysis 

#### 6.1 - Over-representation Analysis (ORA): GO & KEGG 

#### 6.2 - Gene Set Enrichment Analysis (GSEA): GO & KEGG


# Results
### FastQC 

### Salmon Outputs

### Differential Expression Analysis 
#### Wald Test Visualizations 

#### LRT Visualization

### Functional Enrichment Analysis 

#### ORA Visualizations

#### GSEA Visualizations


## Discussion 

## Conclusion 

## References
[1] A. V. Mardanov, M. A. Eldarov, A. V. Beletsky, T. N. Tanashchuk, S. A. Kishkovskaya, and N. V. Ravin, “Transcriptome Profile of Yeast Strain Used for Biological Wine Aging Revealed Dynamic Changes of Gene Expression in Course of Flor Development,” Frontiers in Microbiology, vol. 11, Apr. 2020, doi: https://doi.org/10.3389/fmicb.2020.00538. <br/>
[2] J. Linde et al., “Defining the transcriptomic landscape of Candida glabrata by RNA-Seq,” Nucleic Acids Research, vol. 43, no. 3, pp. 1392–1406, Feb. 2015, doi: https://doi.org/10.1093/nar/gku1357. <br/>
[3] S. Andrews, “FastQC a quality control tool for high throughput sequence data,” Babraham.ac.uk, 2010. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ <br/>
[4] A. Conesa et al., “A survey of best practices for RNA-seq data analysis,” Genome Biology, vol. 17, no. 1, Jan. 2016, doi: https://doi.org/10.1186/s13059-016-0881-8. <br/>
[5] C. Zhang, B. Zhang, L.-L. Lin, and S. Zhao, “Evaluation and comparison of computational tools for RNA-seq isoform quantification,” BMC Genomics, vol. 18, no. 1, Aug. 2017, doi: https://doi.org/10.1186/s12864-017-4002-1. <br/>
[6] R. Patro, G. Duggal, M. I. Love, R. A. Irizarry, and C. Kingsford, “Salmon provides fast and bias-aware quantification of transcript expression,” Nature Methods, vol. 14, no. 4, pp. 417–419, Mar. 2017, doi: https://doi.org/10.1038/nmeth.4197. <br/>
[7] R. Patro, S. M. Mount, and C. Kingsford, “Sailfish enables alignment-free isoform quantification from RNA-seq reads using lightweight algorithms,” Nature Biotechnology, vol. 32, no. 5, pp. 462–464, Apr. 2014, doi: https://doi.org/10.1038/nbt.2862. <br/>
[8] M. I. L. Robinson Charlotte Soneson, Mark D., “Importing transcript abundance with tximport,” bioconductor.org. https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html <br/>
[9] C. Soneson, M. I. Love, and M. D. Robinson, “Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences,” F1000Research, vol. 4, p. 1521, Dec. 2015, doi: https://doi.org/10.12688/f1000research.7563.1. <br/>
[10] M. I. Love, W. Huber, and S. Anders, “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2,” Genome Biology, vol. 15, no. 12, p. 550, Dec. 2014, doi: https://doi.org/10.1186/s13059-014-0550-8. <br/>
[11] M. D. Robinson, D. J. McCarthy, and G. K. Smyth, “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data,” Bioinformatics, vol. 26, no. 1, pp. 139–140, Nov. 2009, doi: https://doi.org/10.1093/bioinformatics/btp616. <br/>
[12] C. W. Law, Y. Chen, W. Shi, and G. K. Smyth, “voom: precision weights unlock linear model analysis tools for RNA-seq read counts,” Genome Biology, vol. 15, no. 2, p. R29, 2014, doi: https://doi.org/10.1186/gb-2014-15-2-r29. <br/>
[13] S. Liu, Z. Wang, R. Zhu, F. Wang, Y. Cheng, and Y. Liu, “Three Differential Expression Analysis Methods for RNA Sequencing: limma, EdgeR, DESeq2,” Journal of Visualized Experiments, no. 175, Sep. 2021, doi: https://doi.org/10.3791/62528.<br/>
[14] G. Yu, L.-G. Wang, Y. Han, and Q.-Y. He, “clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters,” OMICS: A Journal of Integrative Biology, vol. 16, no. 5, pp. 284–287, May 2012, doi: https://doi.org/10.1089/omi.2011.0118. <br/>
[15] A. Subramanian et al., “Gene Set Enrichment analysis: a knowledge-based Approach for Interpreting genome-wide Expression Profiles,” Proceedings of the National Academy of Sciences, vol. 102, no. 43, pp. 15545–15550, Sep. 2005, doi: https://doi.org/10.1073/pnas.0506580102. <br/>
