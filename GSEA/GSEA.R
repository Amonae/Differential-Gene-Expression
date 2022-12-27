library("ggVennDiagram")

# Loading EdgeR and DESeq2 result data

RLE = read.csv("data/EdgeR_RLE_DEGs.csv", row.names = 1)
dim(RLE)
head(RLE)

TMM = read.csv("data/EdgeR_TMM_DEGs.csv", row.names = 1)
dim(RLE)
head(TMM)

UQ = read.csv("data/EdgeR_UQ_DEGs.csv", row.names = 1)
dim(UQ)
head(UQ)

DESeq2 = read.csv("data/DESeq2_DEGs.csv", row.names = 1)
head(DESeq2)

# Isolating DEGs from DESeq2 results (EdgeR results only include DEGs) 

DESeq2 = DESeq2[complete.cases(DESeq2),] # remove NAs 
DESeq2 = DESeq2[DESeq2$padj<0.05,]
summary(DESeq2$log2FoldChange)

#Comparing gene lists using the 3 methods

gene_list = list(row.names(RLE), row.names(TMM), row.names(UQ), row.names(DESeq2))
ggVennDiagram(gene_list, label_alpha = 0, category.names = c('RLE', 'TMM', 'UQ', 'DESeq2'))+ ggplot2::scale_fill_gradient(low="white",high = "green")

# Isolating the 285 genes that are common among all 4

Genes = Reduce(intersect, gene_list)

### Following the protocol on this website https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/


# Installing and loading required packages
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("ggnewscale")
BiocManager::install("ggridges")
library(clusterProfiler)
library(enrichplot)

# Install and load human annotations
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
library(ggplot2)


## Cleaning data. Going to use EdgeR TMM as the base for LFc values
# we want the log2 fold change 
df = TMM[row.names(TMM) %in% Genes,]
gene_list = df$logFC

# name the vector
names(gene_list) = row.names(df)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

### GO Analysis

gse = gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

# Dot plot

require(DOSE)
dotplot(gse)


dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

#GSEA plot
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

# Enrichment map plot

gse_p = pairwise_termsim(gse)
emapplot(gse_p, showCategory = 10)


# Ridgeplot
ridgeplot(gse) + labs(x = "enrichment distribution")

