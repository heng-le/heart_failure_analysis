#Library Installation ####
# installing relevant libaries 
library(tidyverse)
library(ggplot2)
library(ggbiplot)
library(GEOquery)
library(DESeq2)
library(EnhancedVolcano)
library(apeglm)
library(biomaRt)

# Working with metadata ####

# finding metadata on series
data <- getGEO("GSE141910")
data$GSE141910_series_matrix.txt.gz@experimentData

# metadata
clindata <- data[["GSE141910_series_matrix.txt.gz"]]@phenoData@data
t(clindata[1,])

## Let's get the Gene Expression data too
exp <- data[["GSE141910_series_matrix.txt.gz"]]@assayData$exprs

# read in the csv file 
raw_counts <- read.csv("~/Documents/GitHub/heart_failure_analysis/Counts.csv", row.names = 1)

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


#Exploratory Data Analysis ####
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
clindata_hcm <- clindata[clindata$etiology %in% c("Non-Failing Donor","Hypertrophic cardiomyopathy (HCM)"),]
clindata_hcm$etiology <- factor(clindata_hcm$etiology, levels = c("Non-Failing Donor", "Hypertrophic cardiomyopathy (HCM)"))
# then, subset the raw_counts based on the clindata
raw_counts_hcm <- raw_counts[,colnames(raw_counts) %in% rownames(clindata_hcm)]

# now, perform pca on this data
dds_hcm <- DESeqDataSetFromMatrix(countData = raw_counts_hcm,
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
ggplot(df_pc_hcm, aes(x = PC3, y = PC4, label = row.names(df_pc_hcm),color = etiology))+
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
                                 design = ~ etiology + sex)

# we added an interaction term since we want to test whether the effect of sex varies significantly between the two heart conditions 


ddsbat <- DESeq(ddsbat)
resultsNames(ddsbat)

resbat <- results(ddsbat)
resbat_et <- results(ddsbat, name = "etiology_Hypertrophic.cardiomyopathy..HCM._vs_Non.Failing.Donor")
resbat_sex <- results(ddsbat, name = "sex_Male_vs_Female")
#assume no interaction
#combining the two resbats will add up log fold changes 


summary(resbat_et)
summary(resbat_sex)


resbat_et <- resbat_et[order(resbat_et$padj),]
resbat_et
summary(resbat_et) 

# removing all rows with NA 
resbat_et <- na.omit(resbat_et)

# finding the top up/downregulated genes
# For statistically significant upregulated genes
top_upregulated_significant <- resbat_et[resbat_et$padj <= 0.05 & resbat_et$log2FoldChange > 1,]
top_upregulated_significant <- top_upregulated_significant[order(top_upregulated_significant$log2FoldChange, decreasing = TRUE), ][1:5,]

# For statistically significant downregulated genes
top_downregulated_significant <- resbat_et[resbat_et$padj <= 0.05 & resbat_et$log2FoldChange < -1, ]
top_downregulated_significant <- top_downregulated_significant[order(top_downregulated_significant$log2FoldChange, decreasing = FALSE), ][1:5,]

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

resbat_et

# upregulated: BIRC7, SLC1A7, PLD4, AGTR2, TMEM229A
# downregulated: IL1RL1, SAA2, VWA5B1, LCN15, UMODL1-AS1

EnhancedVolcano(resbat_et, x = "log2FoldChange", y="pvalue",  #unshrunken lfc
                lab = rownames(resbat_et),    
                FCcutoff = 1, pCutoff = 0.05,
                #ylim= c(0,10), xlim=c(-5,5),
                selectLab = c("HJHJH"),
                pointSize = 0.5, col=c('grey30', 'red3', 'red3', 'red3'),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                legendPosition = 'bottom',
                title = "Control vs Hypertrophic Cardiomyopathy (HCM)",
                subtitle = "~ etiology + sex")

## Let's also see what's happening with the pvalues & FDR
hist(resbat_et$pvalue, breaks = 0:20/20, col = "grey50", border = "white") ## So few


#################################
# Now, let's move on to GO analysis 
library(clusterProfiler)
library(org.Hs.eg.db)
# read in the annotations
annotations <- read.csv("~/Documents/GitHub/heart_failure_analysis/annotations.csv")

# making a new column in resbat
resbat_et$genes <- row.names(resbat_et)

# making it a dataframe
resbat_et <- as.data.frame(resbat_et)
res_ids <- inner_join(resbat_et, annotations, by=c("genes"="gene_id"))    

allOE_genes <- as.character(res_ids$genes)
head(allOE_genes)


## Extract significant results
sigOE <- dplyr::filter(res_ids, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
sigOE
sigOE_genes <- as.character(sigOE$genes)
head(sigOE_genes)


## Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE_genes,      #significant gene list
                universe = allOE_genes,  #background gene list
                keyType = "ENSEMBL",     #geneID type
                OrgDb = org.Hs.eg.db,    #GO is organism specific
                ont = "BP",              #GO options: BP/MF/CC/ALL (for all three)
                pAdjustMethod = "BH",    # choosing type of multiple testing adjustment
                qvalueCutoff = 0.05, 
                readable = TRUE)


cluster_summary <- data.frame(ego)

#let's look at some columns of the first GO term returned
t(cluster_summary[1,c(1:4,9)])
t(cluster_summary[2,c(1:4,9)])
t(cluster_summary[3,c(1:4,9)])


dotplot(ego, showCategory=20)

OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)

###------------------------------------------------------------------------------

### Functional Class Scoring / Gene Set Enrichment Analysis (GSEA)

## You need then entire gene list analysis with a ranking statistic (preferably logFC, or pvalue)
## Usually run on pathway databases like KEGG/Reactome. But can be run on any gene set database.

## Let's run GSEA using GO BP to see if it looks different

# create a gene list using log2FC - it needs to be a named vector
res_ids <- dplyr::filter(res_ids, genes != "NA") ## remove any NAs
fc_ensembl <- res_ids$log2FoldChange ### Extract the foldchanges
names(fc_ensembl) <- res_ids$genes ### Name each fold change with the corresponding GeneID
fc_ensembl <- sort(fc_ensembl, decreasing = TRUE) ### Sort fold changes in decreasing order
head(fc_ensembl)

gseGO <- gseGO(geneList      = fc_ensembl,
               OrgDb        = org.Hs.eg.db,
               ont          = "BP",
               keyType      = "ENSEMBL",
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

gseGO_res <- gseGO@result ## Extract the GSEA results
head(gseGO_res)
#1 
gseaplot(gseGO, geneSetID = 'GO:0002250')
#2
gseaplot(gseGO, geneSetID = 'GO:0007059')
#conclusion: for patients with HCM, genes related to the modulation of chemical synaptic transmission are upregulated

# You can explore if & how the ORA & GSEA results differ

#### GSEA using KEGG

# KEGG requires Entrez IDs. Let's get the gene list for it

## Remove any NA values (reduces the data by quite a bit)
res_entrez <- dplyr::filter(res_ids, entrezid != "NA")
## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]
foldchanges <- res_entrez$log2FoldChange ### Extract the foldchanges
names(foldchanges) <- res_entrez$entrezid ### Name each fold change
foldchanges <- sort(foldchanges, decreasing = TRUE) ### Sort fold changes in decreasing order
head(foldchanges)

## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms in KEGG
                    minGSSize = 120,  # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result

## Write GSEA results to file
View(gseaKEGG_results)
write.csv(gseaKEGG_results, "gsea_kegg.csv", quote=F)

## Visualizing the GSEA results
## Plot the GSEA plot for a single enriched pathway, `hsa03040`
gseaplot(gseaKEGG, geneSetID = 'hsa04514')

# QUESTION: Is this a positively or negatively enriched pathway?

## We want to look into details of the pathway
# Use the Pathview R package to integrate the KEGG pathway data from clusterProfiler into pathway images

detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts

library(pathview)
## Output images for a single significant KEGG pathway
pathview(gene.data = foldchanges,
         pathway.id = "hsa04514",
         species = "hsa")

## QUESTION: What conclusions can you make from this image?

# Lets try visualize another pathway
pathview(gene.data = foldchanges,
         pathway.id = "hsa04080",
         species = "hsa",
         limit = list(gene = 2, # value gives the max/min limit for foldchanges
                      cpd = 1))

#ClIP-Seq Analysis ####
# RBP - SRSF1 Analysis
srsf1 <- read.csv("~/Documents/GitHub/heart_failure_analysis/srsf1_binding_rbp.csv", header = T)
# plot distribution of binding site records
hist(srsf1$Binding.site.records, breaks = 200, col = "grey50", border = "white", xlim = c(0, 1000))
summary(srsf1$Binding.site.records)
# finding out the 95th and 99th percentile
quantile(srsf1$Binding.site.records, c(0.95, 0.99))
# filter out the top 5% of binding sites
srsf1_top5 <- srsf1[srsf1$Binding.site.records > 273,]
hist(srsf1_top5$Binding.site.records, breaks = 200, col = "grey50", border = "white", xlim = c(0, 1000))

# find out if any of the genes that are differentially expressed overlaps here
srsf1_genes <- srsf1_top5$Target.gene.ID

# combine upregulated and downregulated genes
all_genes <- c(rownames(top_upregulated_significant), rownames(top_downregulated_significant))
# overlap
overlap <- intersect(all_genes, srsf1_genes)

# convert overlap from ensembl to gene names
overlap_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                       filters = "ensembl_gene_id",
                       values = overlap,
                       mart = ensembl)
overlap_genes

# perform gene ontology
ego_overlap <- enrichGO(gene = overlap_genes$ensembl_gene_id,      #significant gene list
                        universe = allOE_genes,  #background gene list
                        keyType = "ENSEMBL",     #geneID type
                        OrgDb = org.Hs.eg.db,    #GO is organism specific
                        ont = "BP",              #GO options: BP/MF/CC/ALL (for all three)
                        pAdjustMethod = "BH",    # choosing type of multiple testing adjustment
                        qvalueCutoff = 0.05, 
                        readable = TRUE)

cluster_summary_clip <- data.frame(ego_overlap)

#let's look at some columns of the first GO term returned
t(cluster_summary_clip[1,c(1:4,9)])
t(cluster_summary_clip[2,c(1:4,9)])
t(cluster_summary_clip[3,c(1:4,9)])


dotplot(ego_overlap, showCategory=20)

OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)