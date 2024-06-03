## ----message=FALSE, warning = FALSE-----------------------------------------------------------------------
# Load necessary libraries
library(survminer)
library(edgeR)
library(EnhancedVolcano)
library(RColorBrewer) # For a colourful plot
library(ggrepel)


## ---------------------------------------------------------------------------------------------------------
# Obtain cutpoint value
COAD.cut <- surv_cutpoint(
  COAD.surv,
  time = "overall_survival",
  event = "event",
  variables = "RAC1")

# Display summary and statistics
summary(COAD.cut)
mean(COAD.surv$RAC1, na.rm = TRUE) 
sd(COAD.surv$RAC1, na.rm = TRUE)
median(COAD.surv$RAC1, na.rm = TRUE)


# Plot cutpoint
plot(COAD.cut, "RAC1", palette = "npg")

# Extract cutpoint
cutpoint <- COAD.cut[["cutpoint"]]$cutpoint


## ---------------------------------------------------------------------------------------------------------
#RAC1 is row 7245
# Assign groups based on cutpoint
groups <- ifelse(counts_matrix_normalised[7245, ] <= cutpoint, "low", "high")
table(groups)


## ---------------------------------------------------------------------------------------------------------
# Normalize data
# We first create a DGEList object
y <- DGEList(counts=counts_matrix, group = groups)
# Setting "low" expression level as control
y$samples$group <- relevel(y$samples$group, ref="low")
# Filter out lowly expressed genes
keep <- filterByExpr(y) 
summary(keep)
y<- y[keep,,keep.lib.sizes=FALSE]


## ----fig.width=50, fig.height=15, warning=FALSE, message=FALSE--------------------------------------------
# Plot unnormalized data
nsamples <- ncol(y) #Number of samples
col <- brewer.pal(nsamples, "Paired")
lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="Log-cpm")

# Performs the trimmed mean of M-values (TMM) normalization calculating normalization factors
y <- calcNormFactors(y, method = "TMM") 

# Creating a design matrix 
design <- model.matrix(~ 0 + y$samples$group)
colnames(design) <- levels(y$samples$group)

#Plotting data after normalisation 
lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Normalised data",ylab="Log-cpm")

# Calculates dispersion using classic mode
y <- estimateDisp(y, design = design, robust = T) # We use robust to treat possible outliers

# Exact test for differential expression
et <- exactTest(y) # Calculates p-values

# Makes decisions about which genes show significant differences
is.de <- decideTestsDGE(et) # Uses p-values
summary(is.de)

# Obtaining a summary of the most significant results in the analysis
# n = Inf => obtaing all results from the analysis
et_tags<-topTags(et, n = Inf)

x<-as.data.frame(et_tags)


# Classify genes into upregulated/stable/downregulated using FDR < 0.05
# logFC >= 1
x$Expression = ifelse(x$FDR < 0.05 & abs(x$logFC) >= 1, 
                      ifelse(x$logFC> 1 ,'Up','Down'),
                      'Stable')
table(x$Expression)

# logFC >= 2
x$Expression2 = ifelse(x$FDR < 0.05 & abs(x$logFC) >= 2, 
                       ifelse(x$logFC> 2 ,'Up','Down'),
                       'Stable')
table(x$Expression2)

# Obtain names of up and down genes
genes_up <- row.names(x[(x$Expression == "Up")&(x$PValue<0.05), ])
genes_down<- row.names(x[(x$Expression == "Down")&(x$PValue<0.05), ])


## ----fig.width=15, fig.height= 10-------------------------------------------------------------------------
# Volcano plot for differential expression
EnhancedVolcano(x,
                lab = rownames(x),
                x= "logFC",
                y= "PValue",
                title = "HIGH Rac vs. LOW Rac",
                pCutoff = 10e-32,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0,
                xlim = c(-5,5),
                # boxedLabels = TRUE,
                col = c("azure4", "palegreen3", "skyblue3", "palevioletred"),
                colAlpha = 0.5,
                legendPosition = 'right')


## ---------------------------------------------------------------------------------------------------------
# Linear modeling for differential expression
v <- voom(y, design, plot=TRUE) #Produces EList object

fit <- lmFit(v, design)
head(coef(fit))

contr <- makeContrasts(high - low, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp) #Multiple testing errors => removes false positives 
top.table <- topTable(tmp, sort.by = "P", n = Inf) #Extract all DEGs, sort by pvalue
head(top.table, 20)

#Get the number of DEGs
length(which(top.table$adj.P.Val < 0.05))
# Classify genes based on logFC and adjust p-values
top.table$Expression = ifelse(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) >= 1, 
                            ifelse(top.table$logFC> 1 ,'Up','Down'),
                            'Stable')
table(top.table$Expression)


## ----fig.width=15, fig.height= 10-------------------------------------------------------------------------
# Volcano plot for differential expression
EnhancedVolcano(top.table,
                lab = top.table$ID,
                x= "logFC",
                y= "adj.P.Val",
                title = "HIGH Rac vs. LOW Rac",
                pCutoff = 10e-32,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0,
                xlim = c(-5,5),
                # boxedLabels = TRUE,
                col = c("azure4", "palegreen3", "skyblue3", "palevioletred"),
                colAlpha = 0.5,
                legendPosition = 'right')

