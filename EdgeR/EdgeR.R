library("edgeR")
library("statmod")


GSE_data = read.csv("data/GSE147507_subset_clean.csv", row.names = 1)
coldata = read.csv("data/design_data.csv", row.names = 1)

dim(GSE_data)
head(GSE_data[,1:8])
head(coldata)

## Creating a DGEList object
# First need to define groups for the samples

groups = factor(paste(coldata$condition,coldata$cell_type,sep="_"), levels = unique(paste(coldata$condition,coldata$cell_type,sep="_")))
coldata$Group = groups # adding "groups" as a column to coldata
head(coldata)

#Creating DGEList
dge = DGEList(GSE_data, group = groups)
dim(dge)
head(dge$samples)


#Going to filter out low count genes 
keep = filterByExpr(dge,group = groups)
dge = dge[keep,keep.lib.size=FALSE]
dim(dge)
head(dge$samples)

# Comparing Different normalization methods
## TMM

dge_tmm = calcNormFactors(dge, method="TMM")
head(dge_tmm$samples)


##RLE

dge_rle = calcNormFactors(dge, method="RLE")
head(dge_rle$samples)


##Upper Quartile
dge_uq = calcNormFactors(dge, method="upperquartile")
head(dge_uq$samples)


#Viewing MDS Plots based on each normalization technique

##TMM
colors = c(rep(c("red", "hotpink", "green4", "green", "blue4", "blue", "purple", "plum"), each = 3), rep(c("gold4", "gold"), each = 2))

plotMDS(dge_tmm, col=colors, pch=16)

legend("topleft", legend= unique(groups), pch=16, col=unique(colors), cex = 0.6)  

##RLE
plotMDS(dge_rle, col=colors, pch=16)

legend("topleft", legend= unique(groups), pch=16, col=unique(colors), cex = 0.6)

##UQ
plotMDS(dge_uq, col=colors, pch=16)

legend("topleft", legend= unique(groups), pch=16, col=unique(colors), cex = 0.6)


# Subsetting and looking at only NHBE cells

NHBE = GSE_data[,1:6]
groups_NHBE = groups[1:6, drop = T]
dge = DGEList(NHBE, group = groups_NHBE)
dim(dge)

## Filtering low count genes
keep = filterByExpr(dge,group = groups_NHBE)
dge = dge[keep,keep.lib.size=FALSE]
dim(dge)
head(dge$samples)


## Calc norm factors
dge_tmm = calcNormFactors(dge, method="TMM")
dge_rle = calcNormFactors(dge, method="RLE")
dge_uq = calcNormFactors(dge, method="upperquartile")


### Estimating dispersion on dge_TMM to start
# edgeR uses the quantile-adjusted conditional maximum likelihood (qCML) method for experiments with single factor

#Dispersion is a measure of the degree of inter-library variation for a gene
# Larger dispersions = higher variance between replicates, which reduce power to detect DE.

# Creating a design matrix 
design = model.matrix(~ 0 + dge_tmm$samples$group)
colnames(design) = levels(dge_tmm$samples$group)

# estimating dispersion
dge_tmm.disp = estimateDisp(dge_tmm, design, robust = T)

# Viewing BCV plot
plotBCV(dge_tmm.disp) # Note there is a downward trend in coefficient of variation as cpm increases, so using common dispersion may not be appropriate


#GLM dispersion estimates

#bin.spline method
dge_tmm_c.disp = estimateGLMCommonDisp(dge_tmm,design, verbose=T) # estimating common dispersion based on glm model
dge_tmm_tr.disp = estimateGLMTrendedDisp(dge_tmm_c.disp,design, method="bin.spline") #trended dispersion
dge_tmm_tag.disp = estimateGLMTagwiseDisp(dge_tmm_tr.disp,design) #tagwise dispersion
plotBCV(dge_tmm_tag.disp)

#DGE Using different dispersion methods
 
# dge_tmm.disp = estimateDisp(dge_tmm, design, verbose=T, robust = T)
fit = glmQLFit(dge_tmm.disp, design)

#Viewing pairwise NHBE Mock vs Cov
NHBE_dge = glmQLFTest(fit, contrast=c(1,-1))
topTags(NHBE_dge, n=10)
NHBE_dge_list = decideTestsDGE(NHBE_dge, adjust.method="BH", p.value=0.05)
EstDisp = row.names(NHBE_dge_list[NHBE_dge_list!= 0,1,1])
sum(abs(NHBE_dge_list))


# dge_tmm_tr.disp = estimateGLMTrendedDisp(dge_tmm_c.disp,design, method="bin.spline")
fit = glmQLFit(dge_tmm_tr.disp, design)

#Viewing pairwise NHBE Mock vs Cov
NHBE_dge = glmQLFTest(fit, contrast=c(1,-1))
topTags(NHBE_dge, n=10)
NHBE_dge_list = decideTestsDGE(NHBE_dge, adjust.method="BH", p.value=0.05)
GLM_Trend = row.names(NHBE_dge_list[NHBE_dge_list!= 0,1,1])
sum(abs(NHBE_dge_list))

# dge_tmm_tag.disp = estimateGLMTagwiseDisp(dge_tmm_tr.disp,design)
fit = glmQLFit(dge_tmm_tag.disp, design)

#Viewing pairwise NHBE Mock vs Cov
NHBE_dge = glmQLFTest(fit, contrast=c(1,-1))
topTags(NHBE_dge, n=10)
NHBE_dge_list = decideTestsDGE(NHBE_dge, adjust.method="BH", p.value=0.05)
GLM_Tag = row.names(NHBE_dge_list[NHBE_dge_list!= 0,1,1])
sum(abs(NHBE_dge_list))

#library("ggVennDiagram")

gene_list = list(EstDisp, GLM_Trend, GLM_Tag)
ggVennDiagram(gene_list, label_alpha = 0, category.names = c('EstDisp', 'GLM_Trend', 'GLM_Tag'))+ ggplot2::scale_fill_gradient(low="white",high = "blue")

### No real difference in dispersion estimate method. Conintuing with EstDisp

# dge_tmm.disp = estimateDisp(dge_tmm, design, verbose=T, robust = T)
#fit = glmQLFit(dge_tmm.disp, design)
#NHBE_dge = glmQLFTest(fit, contrast=c(1,-1))
#topTags(NHBE_dge, n=10)
#NHBE_dge_list = decideTestsDGE(NHBE_dge, adjust.method="BH", p.value=0.05)


#TMM = row.names(NHBE_dge_list[NHBE_dge_list!= 0,1,1])  * Please note, this is the same as the EstDisp object

# Viewing Smear plot


plotSmear(NHBE_dge, de.tags=TMM)
abline(h = c(-2, 2), col = "blue")

# Exploring how normalization affects DGE list in NHBE cells only
## TMM was done above

# RLE
dge1 = estimateDisp(dge_rle,design, robust = T) 
fit = glmQLFit(dge1, design)
NHBE_dge = glmQLFTest(fit, contrast=c(1,-1))
NHBE_dge_list = decideTestsDGE(NHBE_dge, adjust.method="BH", p.value=0.05)
RLE = rownames(dge1)[as.logical(NHBE_dge_list)]
length(RLE) # number of DEGs

plotSmear(NHBE_dge, de.tags=RLE)
abline(h = c(-2, 2), col = "blue")



#UQ
dge1 = estimateGLMCommonDisp(dge_uq,design, verbose=T) # estimating common dispersion based on glm model
dge1 = estimateGLMTrendedDisp(dge1,design, method="bin.spline") #trended dispersion
fit = glmQLFit(dge1, design)
NHBE_dge = glmQLFTest(fit, contrast=c(1,-1))
NHBE_dge_list = decideTestsDGE(NHBE_dge, adjust.method="BH", p.value=0.05)
UQ = rownames(dge1)[as.logical(NHBE_dge_list)]
length(UQ) # number of DEGs

plotSmear(NHBE_dge, de.tags=UQ)
abline(h = c(-2, 2), col = "blue")


# Venn Diagram 

gene_list = list(TMM, RLE, UQ)
ggVennDiagram(gene_list, label_alpha = 0)+ ggplot2::scale_fill_gradient(low="white",high = "blue")

