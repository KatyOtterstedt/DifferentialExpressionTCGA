# DifferentialExpressionTCGA
This repository contains scripts to perform differential expression analysis (DE) and evaluate and summarize the results. It downloads RNA-Seq data from TCGA database and divides patients according to their expression of our gene of interest. The cutpoint to determine this division is obtained through 'Surv_cutpoint' function. 
Differential expression is performed through packages EdgeR and Limma.

Dependencies
------------

- All code is written in R.
- R libraries required include:
  - edgeR
  - limma
  - dplyr
  - survival
  - survminer

Content
------------
| Filename      | Description   | 
|------------------|-------------| 
|1)Preprocessing_DE| This file contains preprocessing of samples to filter samples with low expression, detection of outliers and data normalisation. Should NOT be overlooked. | 
| 2)DE_Analysis      | This file contains the actual differential expression analysis, from obtaining of the cutpoint value to the plotting of differentially expressed genes.| 
