# DifferentialExpressionTCGA
This repository contains scripts and tools for performing differential expression (DE) analysis on TCGA (The Cancer Genome Atlas) data and evaluate and summarize the results. It downloads RNA-Seq data from TCGA database and divides patients according to their expression of our gene of interest. The cutpoint to determine this division is obtained through 'Surv_cutpoint' function. 
Differential expression is performed through packages EdgeR and Limma.

## Table of Contents
* Overview
* Features
* Installation
* Usage

### Overview
The goal of this project is to analyze differential gene expression in cancer patients using TCGA data. The analysis is primarily done using R, leveraging various bioinformatics packages to handle and visualize transcriptomic data.

## Features

- Processing and cleaning TCGA transcriptomic data
- Differential expression analysis
- Data visualization with plots and graphs
- Automation of analysis workflows using R scripts

## Installation and dependencies
To run the scripts in this repository, you need to have R installed on your system along with the required packages. You can install the necessary packages by running:

```R
install.packages(c("tidyverse", "DESeq2", "edgeR", "ggplot2", "limma", "dplyr", "survival", "survminer"))
```

- All code is written in R.
- R libraries required include:
  - edgeR
  - limma
  - dplyr
  - survival
  - survminer

Content
------------
Filtering workflow contains the workflow used when filtering patients before performing the analysis and standard workflow contains the workflow used when filtering patients after performing the analysis. No significant differences were found among these two workflows for the performed analysis, but it cannot be ruled out that differences might arise when performing this analysis using a different cohort or gene of interest.
| Filename      | Description   | 
|------------------|-------------| 
| 1)Preprocessing_DE | This file contains preprocessing of samples to filter samples with low expression, detection of outliers and data normalisation. Should NOT be overlooked. | 
| 2)DE_Analysis      | This file contains the actual differential expression analysis, from obtaining of the cutpoint value to the plotting of differentially expressed genes.| 
