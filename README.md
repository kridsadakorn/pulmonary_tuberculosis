# Correlation between clinical profiles and transmissibility of pulmonary tuberculosis with human and bacterial genetic profiles

## Overview
This repository contains the scripts and codes used for the analysis of whole genome sequencing (WGS) data and phylogenetic analysis as described in our manuscript. The analyses include sequencing, quality filtering, read mapping, variant calling, phylogenetic tree construction, bacterial genetic clustering, host genetic cluster analysis, and statistical analysis.

## Methods
### 1. Whole Genome Sequencing (WGS)
- **DNA Extraction and Library Preparation**: 
  - Genomic DNA was extracted and prepared using the Nextera XT DNA Library Preparation Kit.
  - The library was denatured, diluted, spiked with PhiX control, and sequenced on the Illumina NextSeq 500 platform.
- **Quality Filtering**:
  - Sequencing reads were filtered using Trimmomatic v0.39.
  - Reads shorter than 50 bp were discarded.
  
```
sh 1_1_quality_filtering.sh
```

- **Read Mapping**:
  - Trimmed reads were mapped to the H37Rv reference genome using bwa mem v0.7.17.
  - Duplicate reads were marked using Picard v2.22.8.
  
```
sh 1_2_read_mapping.sh
```

- **Variant Calling**:
  - Variants were called using GATK HaplotypeCaller with a haploid model.
  - Quality control filters were applied to exclude low-quality SNPs and regions.
  
```
sh 1_3_variant_calling.sh
```

### 2. Phylogenetic Analysis
- **Phylogenetic Tree Construction**:
  - A phylogenetic tree was inferred using IQ-TREE v2.0.5 with the GTR+G4 nucleotide substitution model.
  - M. bovis was used as an outgroup to root the tree.
  
```
sh 2_phylogenetic_tree.sh
```

### 3. Bacterial Genetic Clustering
- **Cluster Analysis**:
  - Pairwise SNV differences were counted using MEGA7.
  - Clusters were defined using a 12 SNV threshold and visualized using PopArt.
  
```
sh 3_cluster_analysis.sh
```

### 4. Host Genetic Cluster Analysis
- **SNP Genotyping and Quality Control**:
  - SNP genotyping was performed using the Illumina Human610-Quad BeadChip array.
  - Quality control steps were applied using PLINK.
  
```
sh 4_1_genotype_qc.sh
```

- **Population Stratification**:
  - Principal component analysis and genetic admixture analysis were performed using IPCAPS and ADMIXTURE.
  
```
4_2_pca_admixture.sh
```

  - Create plots for ADMIXTURE and calcuate Fst values for clustering results
  
```
Rscript 4_3_admixture_plot.R
```

### 5. NAT2 Genotyping
- **NAT2 Haplotype-Specific PCR**:
  - HS-PCR was used to identify SNPs in the NAT2 gene region.
  
```
sh 5_nat2_hspcr.sh
```

### 6. Statistical Analysis
- **Statistical Tests**:
  - Statistical analysis was performed using PASW Statistics 21.0.
  - Various tests, including chi-square, Fisher’s exact test, and logistic regression, were used to analyze the data.

```
Rscript 6_statistical_analysis.R
```

## Files
- `1_1_quality_filtering.sh`: Script for quality filtering of sequencing reads using Trimmomatic.
- `1_2_read_mapping.sh`: Script for mapping reads to the reference genome using bwa mem.
- `1_3_variant_calling.sh`: Script for variant calling using GATK HaplotypeCaller.
- `2_phylogenetic_tree.sh`: Script for constructing the phylogenetic tree using IQ-TREE.
- `3_cluster_analysis.sh`: Script for bacterial genetic clustering using MEGA7 and PopArt.
- `4_1_genotype_qc.sh`: Script for quality control of SNP genotyping data using PLINK.
- `4_2_pca_admixture.sh`: Script for population stratification analysis using IPCAPS and ADMIXTURE.
- `4_3_admixture_plot.R`: Script for creating plots for the results of ADMIXTURE and calcuate Fst values for clustering results.
- `5_nat2_hspcr.sh`: Script for NAT2 genotyping using HS-PCR.
- `6_statistical_analysis.R`: R script for statistical analysis using PASW Statistics.

## Contact
For any questions or issues, please contact the corresponding author of the manuscript.
