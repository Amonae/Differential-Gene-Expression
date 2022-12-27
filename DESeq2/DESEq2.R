library(DESeq2)
BiocManager::install("apeglm")

head(GSEA_data[,1:8])
head(coldata)
coldata$condition = factor(coldata$condition, levels = c("mock", "covid"))

# Doing DGE analysis for only NHBE cells 

dds_NHBE = DESeqDataSetFromMatrix(countData = GSEA_data[,1:6],
                                  colData = coldata[1:6,],
                                  design = ~ condition)    # Creating DESeq dataset


dds_NHBE = DESeq(dds_NHBE)   # Calling DESeq function

#Results object
dge_NHBE = results(dds_NHBE, alpha = 0.05, lfcThreshold = 0)
dge_NHBE    # Note that the LFC is based on covid vs mock. So positive values represent genes that are upregulated in cells with covid


ordered_p1 = dge_NHBE[order(dge_NHBE$pvalue),]
ordered_p1

# Going to shrink for better visualization and ranking of genes.

dge_NHBE_shrink = lfcShrink(dds_NHBE, coef="condition_covid_vs_mock", type="apeglm")
dge_NHBE_shrink

mcols(dge_NHBE_shrink)$description

summary(dge_NHBE_shrink)  # note that here, the minimum pvalue is < 0.1

summary(dge_NHBE_shrink, alpha = 0.05) # setting the minimum pvalue to < 0.05

sum(dge_NHBE_shrink$padj < 0.05, na.rm=TRUE) # reports how many genes have padj < 0.05

ordered_p2 = dge_NHBE_shrink[order(dge_NHBE_shrink$pvalue),] # This orders the output by lowest pvalues

head(ordered_p2, 7) 

plotMA(dge_NHBE_shrink, ylim=c(-2,2))  # This plot shows LFC of genes with significant genes colored blue.


# data needs to be transformed for PCA to remove mathematically the sources of unwanted variations
rld = rlog(dds, blind=FALSE)
plotPCA(rld)

## Saving results to csv
DESeq2 = data.frame(dge_NHBE_shrink)
write.csv(DESeq2, file = "data/DESeq2_DEGs.csv")

