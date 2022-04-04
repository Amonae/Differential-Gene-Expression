library("edgeR")
library("statmod")

#Section2

# For an initial look at the data check how many reads there are in each sample (column) and present this as a bar graph (Figure 1). Explain why some samples have more reads than others and the 

raw_counts_df = read.csv(file.choose(), header = T, row.names = 1)
View(raw_counts_df)

totals = colSums(raw_counts_df) #### sum counts for each column
colors <- c(rep(c("red", "hotpink", "green", "green4", "blue", "blue4", "purple", "plum"), each = 3), rep(c("gold", "gold4"), each = 2))
barplot(totals, col = colors)
# Load the count matrix into EdgeR and perform normalization using the TMM, RLE and upper quartile methods

groups = c(rep(c("Series 1 Mock", "Series 1 SARS", "Series 2 Mock", "Series 2 SARS", "Series 6 Mock", "Series 6 SARS", "Series 7 Mock", "Series 7 SARS"), each = 3), rep(c("Series 15 Healthy", "Series 15 Covid"), each = 2))
#groups
dge = DGEList(raw_counts_df, group = groups)  ### create DGElist
#dge

tmm_normalization = calcNormFactors(dge, method="TMM")  #TMM normalization factors. this is a dgelist
#tmm_normalization
#pseudo_tmm <- log2(cpm(tmm_normalization) + 1) # pseudocounts tmm. this is a matrix
#View(pseudo_tmm)

rle_normalization = calcNormFactors(dge, method="RLE")  #RLE normalization factors. this is a dgelist
#rle_normalization
#pseudo_rle <- log2(cpm(rle_normalization) + 1) # pseudocounts rle. this is a matrix
#pseudo_rle

uq_normalization = calcNormFactors(dge, method="upperquartile")  #Upper Quartile normalization factors. this is a dgelist
#uq_normalization
#pseudo_uq <- log2(cpm(uq_normalization) + 1) # pseudocounts uq. this matrix
#pseudo_uq


# Provide guidelines for how to generate a PCA plot that clearly depicts the relationships between the different samples based on the TMM normalised count matrix.

colors <- c(rep(c("red", "hotpink", "green", "green4", "blue", "blue4", "purple", "plum"), each = 3), rep(c("gold", "gold4"), each = 2))

plotMDS(tmm_normalization, col=colors, pch=16)

legend("topleft", legend=c("Series 1 Mock", "Series 1 SARS", "Series 2 Mock", "Series 2 SARS", "Series 6 Mock", "Series 6 SARS", "Series 7 Mock", "Series 7 SARS", "Series 15 Healthy", "Series 15 Covid"), pch=16, col=c("red", "hotpink", "green", "green4", "blue", "blue4", "purple", "plum", "gold", "gold4"), cex = 0.6)  ##Make legend separate from plot






#Section 3

# Repeat the analysis of differential gene expression, but just between the mock and SARS-CoV-2 infected NHBE samples that you analysed in section 2 above. 


groups_NHBE = as.factor(rep(c("NHBE.Mock", "NHBE.SARS"), each = 3)) # define groups
#str(groups_NHBE)

design <- model.matrix(~0+groups_NHBE)   #making a design matrix. This design matrix simply links each group to the samples that belong to it. Each row of the design matrix corresponds to a sample whereas each column represents a coefficient corresponding to one of the 2 groups.

colnames(design) <- levels(groups_NHBE)
#design

NHBE_raw = raw_counts_df[c(1:6)]        #create NHBE dataframe
#View(NHBE_raw)

NHBE_dge = DGEList(NHBE_raw, group = groups_NHBE)  ### create DGElist
#View(NHBE_dge)

NHBE_tmm = calcNormFactors(NHBE_dge, method="TMM")  #TMM normalization factors. this is a dgelist
NHBE_tmm = estimateDisp(NHBE_tmm, design, robust=TRUE)  ## This returns a DGEList object with additional components (common.dispersion, trended.dispersion and tagwise.dispersion) added to hold the estimated dispersions. 

NHBE_rle = calcNormFactors(NHBE_dge, method="RLE")  #RLE normalization factors. this is a dgelist
NHBE_rle = estimateDisp(NHBE_rle, design, robust=TRUE)

NHBE_uq = calcNormFactors(NHBE_dge, method="upperquartile")  #UQ normalization factors. this is a dgelist
NHBE_uq = estimateDisp(NHBE_uq, design, robust=TRUE)


#Compare the lists of altered genes detected by each method and present the overlaps as a Venn diagram (Figure 3).


NHBE_contrast = makeContrasts(NHBE.Mock-NHBE.SARS, levels=design) #contrast the 2 groups. In subsequent results, a positive log2-fold-change (logFC) will indicate a gene up-regulated Mock relative to SARS cells, whereas a negative logFC will indicate a gene more highly expressed in SARS cells 

NHBE_TMM_fit <- glmQLFit(NHBE_tmm, design, robust=TRUE)  #The estimation of QL dispersions is performed using the glmQLFit function. This returns a DGEGLM object with the estimated values of the GLM coefficients for each gene

DE_TMM = glmQLFTest(NHBE_TMM_fit, contrast=NHBE_contrast) # generates a table with logFC, logCPM, F, Pval, and FDR   TMM
DE_TMM_adj = decideTestsDGE(DE_TMM) ### table says which genes are DEG. 1 = yes 0 = no

#topTags(DE_TMM)
summary(DE_TMM_adj)
#View(DE_TMM_adj)

NHBE_RLE_fit <- glmQLFit(NHBE_rle, design, robust=TRUE)  
DE_RLE = glmQLFTest(NHBE_RLE_fit, contrast=NHBE_contrast) #   RLE
DE_RLE_adj = decideTestsDGE(DE_RLE)
summary(DE_RLE_adj)

NHBE_UQ_fit <- glmQLFit(NHBE_uq, design, robust=TRUE)  
DE_UQ = glmQLFTest(NHBE_UQ_fit, contrast=NHBE_contrast) #   upperquartile
DE_UQ_adj = decideTestsDGE(DE_UQ)
summary(DE_UQ_adj)

comparison = cbind(DE_TMM_adj, DE_RLE_adj, DE_UQ_adj)  #add all results together: TMM RLE UQ

venn_comp = vennCounts(comparison, include = 'both')
venn_diag = vennDiagram(venn_comp, names = c('TMM', "RLE", "UQ"))

#Now use 3 different software packages, EdgeR, DeSeq2 and Limma-Voom, and Log2FC > 1 or <-1, p adjusted-value < 0.05. Use the default parameters (ie use EdgeR with TMM normalization) and note these in your answer

#EdgeR

LFC1_TMM_fit <- glmQLFit(NHBE_tmm, design, robust=TRUE)
DE_TMM_LFC1 = glmQLFTest(LFC1_TMM_fit, contrast=NHBE_contrast) 
DE_TMM_LFC1_adj = decideTestsDGE(DE_TMM_LFC1, lfc = 1) ###               LFC = 1
summary(DE_TMM_LFC1_adj)

LFCn1_TMM_fit <- glmQLFit(NHBE_tmm, design, robust=TRUE)
DE_TMM_LFCn1 = glmQLFTest(LFCn1_TMM_fit, contrast=NHBE_contrast) 
DE_TMM_LFCn1_adj = decideTestsDGE(DE_TMM_LFCn1, lfc = -1) ###               LFC = -1 ***don't think it's necessary to do both
summary(DE_TMM_LFCn1_adj)

#DESeq2

library(DESeq2)
library(DEFormats)

NHBE_dds = as.DESeqDataSet(NHBE_dge)  #convert EdgeR data to DESeq2 object
NHBE_dds_norm = DESeq(NHBE_dds)  # calculate normalization of DESeq2 data

#sizeFactors(NHBE_dds_norm)   # view normalization factors
#results(NHBE_dds_norm)
#dds_norm_to_dge = as.DGEList(NHBE_dds_norm)

NHBE_dds_adj = results(NHBE_dds_norm, contrast = NHBE_contrast, alpha = 0.05, lfcThreshold = 1) ###dataframe with lfc cutoff 1 and pvalue cutoff 0.05

#mcols(NHBE_dds_adj, use.names=T)
#NHBE_dds_adj %>% data.frame() %>% View()
summary(NHBE_dds_adj)

NHBE_dds_adj_df = as.data.frame(NHBE_dds_adj)     # convert DESeq2 results to DF
NHBE_dds_adj_1 = NHBE_dds_adj_df[c(2,6)]        # select pvalue and lfc columns

NHBE_dds_adj_1[is.na(NHBE_dds_adj_1)] = 0      # replace NAs with 0
View(NHBE_dds_adj_1)

for (i in 1:nrow(NHBE_dds_adj_1)) {         #convert DEG LFC values to 1 or 0 ** LFC >1 = 1 else it's 0
  if (NHBE_dds_adj_1[i,1] > 1){
    NHBE_dds_adj_1[i,1] = 1}
  else{
    NHBE_dds_adj_1[i,1] = 0
  }
}

for (i in 1:nrow(NHBE_dds_adj_1)) {         #convert DEG p values to 1 or 0. 
  if (NHBE_dds_adj_1[i,2] < 0.5){
    NHBE_dds_adj_1[i,2] = 1}
  else {NHBE_dds_adj_1[i,2] = 0}
}


add_gene_1 = matrix(nrow = 21797)                # new empty matrix

for (i in 1:nrow(NHBE_dds_adj_1)) {         #add 1 for values that have the right p and lfc value 
  if (NHBE_dds_adj_1[i,1] == 1 & NHBE_dds_adj_1[i,2] == 1){
    add_gene_1[i,1] = 1}
  else {add_gene_1[i,1] = 0}
}
colnames(add_gene_1) = "DESeq.LFC.1"
View(add_gene_1)

NHBE_dds_adj_df = as.data.frame(NHBE_dds_adj)     # convert DESeq2 results to DF
NHBE_dds_adj_n1 = NHBE_dds_adj_df[c(2,6)]        # select pvalue and lfc columns

NHBE_dds_adj_n1[is.na(NHBE_dds_adj_n1)] = 0      # replace NAs with 0
View(NHBE_dds_adj_n1)

for (i in 1:nrow(NHBE_dds_adj_n1)) {         #convert DEG LFC values to 1 or 0 ** both up and downregulated  have value 1
  if (NHBE_dds_adj_n1[i,1] < -1){
    NHBE_dds_adj_n1[i,1] = 1}
  else{
    NHBE_dds_adj_n1[i,1] = 0
  }
}



for (i in 1:nrow(NHBE_dds_adj_n1)) {         #convert DEG p values to 1 or 0. 
  if (NHBE_dds_adj_1[i,2] < 0.5){
    NHBE_dds_adj_n1[i,2] = 1}
  else {NHBE_dds_adj_n1[i,2] = 0}
}


add_gene_n1 = matrix(nrow = 21797)                # new empty matrix

for (i in 1:nrow(NHBE_dds_adj_n1)) {         #add 1 for values that have the right p and lfc value 
  if (NHBE_dds_adj_n1[i,1] == 1 & NHBE_dds_adj_n1[i,2] == 1){
    add_gene_n1[i,1] = 1}
  else {add_gene_n1[i,1] = 0}
}
colnames(add_gene_n1) = "DESeq.LFC.n1"
View(add_gene_n1)


comparison_2 = cbind(DE_TMM_LFC1_adj, add_gene)  #add DESEQ2 and EdgeR results together
View(comparison_2)

### Limma Voom

#design <- model.matrix(~0+groups_NHBE)  ### This is already defined above. Just a reminder

NHBE_tmm_voom = calcNormFactors(NHBE_dge, method="TMM")  #TMM normalization factors. this is a dgelist
#NHBE_tmm_voom = estimateDisp(NHBE_tmm_voom, design, robust=TRUE)  ## This returns a DGEList object with additional components 

NHBE_voom = voom(NHBE_tmm_voom, design = design) # Voom analysis
NHBE_voom_fit <- lmFit(NHBE_voom, design)                              # fit linear model
NHBE_voom_contrasts <- contrasts.fit(NHBE_voom_fit, NHBE_contrast)      #contrasts for each gene
NHBE_voom_adj_1 = decideTests(NHBE_voom_contrasts, lfc = 1) ###               LFC = 1 decidetest
NHBE_voom_adj_n1 = decideTests(NHBE_voom_contrasts, lfc = -1) ###               LFC = -1 decidetest
summary(NHBE_voom_adj_1)
summary(NHBE_voom_adj_n1)




##### Combining everything

comparison_up = cbind(DE_TMM_LFC1_adj, NHBE_voom_adj_1, add_gene_1)#add DESEQ2 and EdgeR results together

comparison_down = cbind(DE_TMM_LFCn1_adj,  NHBE_voom_adj_n1, add_gene_n1)

venn_comp_up = vennCounts(comparison_up, include = 'both')
venn_diag_up = vennDiagram(venn_comp_up)

venn_comp_down = vennCounts(comparison_down, include = 'both')
venn_diag_down = vennDiagram(venn_comp_down)


