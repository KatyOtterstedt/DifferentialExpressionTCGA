## ----message=FALSE, warning = FALSE--------------
# Load necessary libraries
library(RTCGA.clinical)
library(TCGAbiolinks)
library(tidyverse)
library(edgeR)
library(SummarizedExperiment)


## ----message=FALSE, warning = FALSE--------------
# Fetch clinical data for COAD
clinical_COAD <- GDCquery_clinic("TCGA-COAD")
barcodes_to_download <- clinical_COAD$bcr_patient_barcode


## ----message=FALSE, warning = FALSE, eval = FALSE----
# Fetch and prepare gene expression data
query_TCGA <- GDCquery(project = "TCGA-COAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts",
                       sample.type = "Primary Tumor",
                       barcode= barcodes_to_download,
                       access = "open")
 
GDCdownload(query_TCGA)
tcga_coad_SE <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
 
#Extract counts and gene metadata
counts_matrix <- assay(tcga_coad_SE)
gene_metadata <- as.data.frame(rowData(tcga_coad_SE))
row.names(counts_matrix) <- (gene_metadata$gene_name)
 

## ------------------------------------------------
# Filter rows based on gene names (filtering genes that are miRNA or lncRNA
counts_matrix_filt_bef <- counts_matrix[!grepl("^MIR|^LINC|^Y_", rownames(counts_matrix)), ]

# Check dimensions of matrices
dim(counts_matrix)
dim(counts_matrix_filt_bef)



## ----message=FALSE, warning = FALSE--------------
# Filter and normalize data
source("tcga_replicateFilter.R") #Contains function that filters patients according to their barcode
total_patients <- c(colnames(counts_matrix))
filtered_patients <- tcga_replicateFilter(tsb = total_patients,
                                          analyte_target = "RNA")
counts_matrix <- counts_matrix[,filtered_patients]


## ----message=FALSE, warning = FALSE--------------
counts_matrix_normalised <- cpm(counts_matrix)


## ----message=FALSE, warning = FALSE--------------
counts_matrix_normalised <- as.data.frame(counts_matrix_normalised) 

# Shorten patient names for merging
long_patient_names <- colnames(counts_matrix_normalised)

shortened_patient_names <- substr(long_patient_names, 1, 12)

# Converting from vector to df to be added as a row to our normalised counts matrix 
shortened_patient_names <- as.data.frame(t(shortened_patient_names))

colnames(shortened_patient_names) <- colnames(counts_matrix_normalised)

# Adding short names as a row to normalised matrix 
counts_matrix_normalised <- rbind(shortened_patient_names, counts_matrix_normalised)

# Changing the name to make the merge easier
row.names(counts_matrix_normalised)[1] <- "bcr_patient_barcode"

# Check gene index and missing values
which(row.names(counts_matrix_normalised) %in% c("RAC1"))
table(is.na(counts_matrix_normalised[7245,])) 


## ----message=FALSE, warning = FALSE--------------
# Create event variable for survival analysis
clinical_COAD$event <- ifelse(clinical_COAD$vital_status == "Alive",
                                 0,
                                 1)
table(clinical_COAD$event) 

#Create an overal survival variable
clinical_COAD$overall_survival <- ifelse(clinical_COAD$vital_status == "Alive",
                                         clinical_COAD$days_to_last_follow_up, #If alive => last follow up to OS
                                         clinical_COAD$days_to_death) #If dead => death to OS



## ----message=FALSE, warning = FALSE--------------
# Create event variable for 5-year survival analysis
clinical_COAD$event_5y <- ifelse(clinical_COAD$overall_survival >= 1825,
                                 0, 
                                 clinical_COAD$event)
table(clinical_COAD$event_5y)

# Adjust overall survival for 5-year analysis
clinical_COAD$overall_survival_5y <- ifelse(clinical_COAD$overall_survival >= 1825,
                                            1825,
                                            clinical_COAD$overall_survival)


## ----message=FALSE, warning = FALSE--------------
# Select relevant columns for survival analysis
which(colnames(clinical_COAD) %in% c("overall_survival_5y", "event_5y", "bcr_patient_barcode",
                                     "overall_survival", "event"))

COAD.surv <- clinical_COAD[,c(70,71,72,73,74)]
table(is.na(COAD.surv$event_5y))
COAD.surv <- COAD.surv %>%
  filter(!is.na(event_5y))

# Merge survival data with expression data for gene of interest
COAD.surv <- merge(COAD.surv, t(counts_matrix_normalised[c(1,7245),]), by = "bcr_patient_barcode")
# Converting the column to type integer 
COAD.surv$RAC1 <- as.integer(COAD.surv$RAC1)
table(is.na(COAD.surv$RAC1))


## ----message=FALSE, warning = FALSE--------------
# Visualize expression distribution
boxplot(COAD.surv$RAC1, col = "palevioletred")

# Identify and handle outliers
# Get mean and standard deviation
mean = mean(COAD.surv$RAC1)
std = sd(COAD.surv$RAC1)

# Get threshold values for outliers
Tmin = mean-(3*std)
Tmax = mean+(3*std)

# Find outlier expression values
outliers <- which(COAD.surv$RAC1 < Tmin | COAD.surv$RAC1>= Tmax)

# Assing NA values to outliers
COAD.surv$RAC1[outliers] <- NA

# Convert RAC1 column to integer
COAD.surv$RAC1 <- as.integer(COAD.surv$RAC1)

