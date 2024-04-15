# installing relevant libaries 
library(tidyverse)
library(ggplot2)
library(ggbiplot)
library(GEOquery)
library(DESeq2)
library(EnhancedVolcano)
library(apeglm)
library(biomaRt)


# finding metadata on series
data <- getGEO("GSE141910")
data$GSE141910_series_matrix.txt.gz@experimentData

# metadata
clindata <- data[["GSE141910_series_matrix.txt.gz"]]@phenoData@data
t(clindata[1,])

## Let's get the Gene Expression data too
exp <- data[["GSE141910_series_matrix.txt.gz"]]@assayData$exprs

# read in the csv file 
raw_counts <- read.csv("Counts.csv", row.names = 1)

rownames(clindata) <- clindata$title

# making sure that names match 
all(rownames(clindata) %in% colnames(raw_counts))

# picking columns of interest 
clindata <- clindata[,43:46]

# renaming columns 
colnames(clindata)[colnames(clindata) == "age:ch1"] <- "age"
colnames(clindata)[colnames(clindata) == "etiology:ch1"] <- "etiology"
colnames(clindata)[colnames(clindata) == "race:ch1"] <- "race"
colnames(clindata)[colnames(clindata) == "Sex:ch1"] <- "sex"

head(clindata)

# create a new categorical column for age, where < 55 is "young", and >= 55 is "old"
clindata$age_cat <- ifelse(clindata$age < 55, "young", "old")

# make sure that the grouping variables are factors 
clindata$etiology <- as.factor(clindata$etiology)
clindata$race <- as.factor(clindata$race)
clindata$sex <- as.factor(clindata$sex)
clindata$age_cat <- as.factor(clindata$age_cat)

barplot(table(rowSums(raw_counts>0)))
abline(v=83, col = "red", lwd=2) #Half of the non-heart-failure samples
abline(v=ncol(raw_counts)/2, col = "orange", lwd=2) #Half of all samples


keep <- rowSums(raw_counts > 0) >= ncol(raw_counts)/2
table(keep) ## So many rows are false! Can throw out these genes
raw_counts <- raw_counts[keep,]
dim(raw_counts) #Left with 10,976 genes from 35784 starting genes

## Now, let's remove samples with >70% of genes with 0 counts - poor samples
## see how many genes per sample have 0 reads
zero <- colSums(raw_counts == 0)/nrow(raw_counts) #proportion of genes with 0 expression per sample
hist(zero, breaks = 50)
abline(v = 0.7, col="red") #line at samples with 70% gene with 0 counts
# in this case, the maximum proportion of genes with 0 counts is 0.32
summary(zero)

good_samples <- zero <0.7
raw_counts <- raw_counts[,good_samples]


# subset the clindata for control and hcm
clindata_hcm <- clindata[clindata$etiology %in% c("Non-Failing Donor", "Dilated cardiomyopathy (DCM)"),]

# then, subset the raw_counts based on the clindata
raw_counts_hcm <- raw_counts[,colnames(raw_counts) %in% rownames(clindata_hcm)]

# now, perform pca on this data
dds_hcm <- DESeqDataSetFromMatrix(countData = raw_counts_hcm,, 
                              colData = clindata_hcm,    
                              design = ~1)   
vsd_hcm <- vst(dds_hcm) 
vsd.dat_hcm <- as.data.frame(assay(vsd_hcm))
pc_vsd_hcm <- prcomp(t(vsd.dat_hcm)) 
pc.var_hcm <- pc_vsd_hcm$sdev^2
pc.var.per_hcm <- round(pc.var_hcm/sum(pc.var_hcm)*100,1) #also rounding off the percentage to 1 decimal
head(pc.var.per_hcm)

# PCA Analysis scree plot
ggscreeplot(pc_vsd_hcm)+
  geom_col(fill = 'cornflowerblue') + 
  geom_line()+
  theme_bw()

pc.data_hcm <- as.data.frame(pc_vsd_hcm$x)
df_pc_hcm <- cbind(clindata_hcm, pc.data_hcm) # Let's cbind to original df so we have metadata
head(df_pc_hcm)[1:10]

# PCA plot
ggplot(df_pc_hcm, aes(x = PC1, y = PC2, label = row.names(df_pc_hcm),color = etiology))+
  #geom_text()+ ##Uncomment If you wanna check which is the outlier
  geom_point(size=2)+
  xlim(c(-150,150))+ ylim(c(-50,50))+
  theme_bw()

# checking for balance
# clindata_hcm %>%
#   group_by(sex, etiology) %>%
#   dplyr::summarise(frequency = n()) 

ddsbat <- DESeqDataSetFromMatrix(raw_counts_hcm,
                                 colData = clindata_hcm,
                                 design = ~ etiology + sex + etiology:sex)

# we added an interaction term since we want to test whether the effect of sex varies significantly between the two heart conditions 


ddsbat <- DESeq(ddsbat)
resultsNames(ddsbat)

resbat <- results(ddsbat)
resbat <- resbat[order(resbat$padj),]
resbat
summary(resbat) 

# removing all rows with NA 
resbat <- na.omit(resbat)

# finding the top up/downregulated genes
# For statistically significant upregulated genes
top_upregulated_significant <- resbat[resbat$padj <= 0.05 & resbat$log2FoldChange > 1, ]
top_upregulated_significant <- top_upregulated_significant[order(top_upregulated_significant$log2FoldChange, decreasing = TRUE), ][1:5, ]

# For statistically significant downregulated genes
top_downregulated_significant <- resbat[resbat$padj <= 0.05 & resbat$log2FoldChange < -1, ]
top_downregulated_significant <- top_downregulated_significant[order(top_downregulated_significant$log2FoldChange, decreasing = FALSE), ][1:5, ]

# using biomart to find out gene names 
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = 'asia')
downregulated_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
               filters = "ensembl_gene_id",
               values = rownames(top_downregulated_significant),
               mart = ensembl)

upregulated_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                             filters = "ensembl_gene_id",
                             values = rownames(top_upregulated_significant),
                             mart = ensembl)

upregulated_genes
downregulated_genes

resbat

# upregulated: FMO5, ADAMTS16, OVOS2, IGHGP
# downregulated: UPK1B, PRG4, CXCL1, CFB, ANKRD26P1

EnhancedVolcano(resbat, x = "log2FoldChange", y="pvalue",  #unshrunken lfc
                lab = rownames(resbat),
                selectLab = c('FMO5','ADAMTS16','OVOS2','IGHGP','UPK1B',  # top 5up DEG
                              'PRG4','CXCL1','CFB','ANKRD26P1'),          # genes of interest
                FCcutoff = 1, pCutoff = 0.05,
                #ylim= c(0,10), xlim=c(-5,5),
                pointSize = 0.5, col=c('grey30', 'red3', 'red3', 'red3'),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                legendPosition = 'bottom',
                title = "control_vs_dcm",
                subtitle = "~ etiology + sex + etiology:sex")

## Let's also see what's happening with the pvalues & FDR
hist(resbat$pvalue, breaks = 0:20/20, col = "grey50", border = "white") ## So few


#################################
# Now, let's move on to GO analysis 
library(clusterProfiler)
library(org.Hs.eg.db)
# read in the annotations
annotations <- read.csv("annotations.csv")

# making a new column in resbat
resbat$genes <- row.names(resbat)

# making it a dataframe
resbat <- as.data.frame(resbat)
res_ids <- inner_join(resbat, annotations, by=c("genes"="gene_id"))    

allOE_genes <- as.character(res_ids$genes)
head(allOE_genes)


## Extract significant results
sigOE <- dplyr::filter(res_ids, padj < 0.05)
sigOE
sigOE_genes <- as.character(sigOE$genes)
head(sigOE_genes)


## Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE_genes,      #significant gene list
                universe = allOE_genes,  #background gene list
                keyType = "ENSEMBL",     #geneID type
                OrgDb = org.Hs.eg.db,    #GO is organism specific
                ont = "ALL",              #GO options: BP/MF/CC/ALL (for all three)
                pAdjustMethod = "BH",    # choosing type of multiple testing adjustment
                qvalueCutoff = 0.05, 
                readable = TRUE)


cluster_summary <- data.frame(ego)

#let's look at some columns of the first GO term returned
t(cluster_summary[1,c(1:4,9)])
