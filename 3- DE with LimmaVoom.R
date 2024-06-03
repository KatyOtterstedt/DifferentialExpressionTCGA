# Beginning differential expression analysis
library(edgeR)
y <- DGEList(counts=counts_matrix, group = groups)
# Setting "low" expression level as control
y$samples$group <- relevel(y$samples$group, ref="low")
keep <- filterByExpr(y) # Filter out lowly expressed genes
summary(keep)
# Mode   FALSE    TRUE 
# logical   41317   19343
y<- y[keep,,keep.lib.sizes=FALSE]

# Performs the TMM normalization
y<- calcNormFactors(y, method = "TMM") #calcNormFactors doesn’t normalize the data,
#it just calculates normalization factors for use downstream
#the samples for most genes. The default method for computing these scale factors
#uses a trimmed mean of M-values (TMM)

plotMDS(y, labels = NULL, col = as.numeric(y$samples$group))

#---Differential expression
#voom transformation is applied to the normalized and filtered DGEList object
#voom tmb recibe raw counts para normalizar sino
v <- voom(y, design, plot=TRUE) #Produces EList object

fit <- lmFit(v, design)
head(coef(fit))

contr <- makeContrasts(high - low, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp) #Multiple testing errors => removes false positives 
top.table <- topTable(tmp, sort.by = "P", n = Inf) #Extract all DEGs, sort by pvalue
head(top.table, 20)


length(which(top.table$adj.P.Val < 0.05)) #Gets the number of DEGs
#6976
x_et$ID <- row.names(x_et)
common_genes <- intersect(top.table$ID, row.names(x_et))

# Obtener los nombres de las filas en común

common_merged<-merge(top.table, x_et,by= "ID", all = TRUE) # FALSE so that only rows with data from both df are included in the output
common_merged <- common_merged[,c(1,2,8)]
colnames(common_merged)[2] <- "logFC Limma Voom"
colnames(common_merged)[3] <- "logFC EdgeR"


top.table$Expression = ifelse(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) >= 1, 
                            ifelse(top.table$logFC> 1 ,'Up','Down'),
                            'Stable')
table(top.table$Expression)
# Down Stable     Up 
# 21  19326     22 

setwd("C:/Users/byb/Documents/Bioinfo Tesina/Differential Expression")
png("Volcano_Limma_Voom.png", width = 3000, height = 2000, res = 300)
EnhancedVolcano(top.table,
                lab = top.table$ID,
                x= "logFC",
                y= "adj.P.Val",
                title = "HIGH Rac vs. LOW Rac",
                #Los nuestros de interés dan todos  no sig y no FC
                # selectLab = c('VAV1','VAV2','VAV3',
                #               'PREX1','PREX2','TIAM1','TIAM2',
                #               'RASGRF2'),
                pCutoff = 10e-32,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0,
                xlim = c(-5,5),
                # boxedLabels = TRUE,
                col = c("azure4", "palegreen3", "skyblue3", "palevioletred"),
                colAlpha = 0.5,
                legendPosition = 'right')
dev.off()
