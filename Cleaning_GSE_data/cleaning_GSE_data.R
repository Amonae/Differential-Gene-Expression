# GSE data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507 human data
# Downloaded and unzipped file on local computer

GSE_data = read.table(file = file.choose(), sep = '\t', header = TRUE, row.names = 1)  
dim(GSE_data) # 21797   78


Series = paste("Series", c(1,2,6,7,15), "_", sep = "")
GSE_data = GSE_data %>% select(starts_with(Series))


#only want a subset of this data

#Renaming columns for clarity
columns = c(paste("S1_Mock", 1:3, sep = ""), paste("S1_Cov", 1:3, sep = ""), paste("S2_Mock", 1:3, sep = ""),paste("S2_Cov", 1:3, sep = ""), paste("S6_Mock", 1:3, sep = ""),paste("S6_Cov", 1:3, sep = ""), paste("S7_Mock", 1:3, sep = ""),paste("S7_Cov", 1:3, sep = ""), paste("S15_Mock", 1:2, sep = ""),paste("S15_Cov", 1:2, sep = ""))

colnames(GSE_data) = columns 
head(GSE_data)

write.csv(GSE_data, file = "GSE147507_subset_clean.csv")

# Need to make a df outlining the experimental design

coldata = colnames(GSE_data)
condition = c(rep(c("mock", "covid"), each = 3,4),rep(c("mock", "covid"), each = 2))
condition= factor(condition, levels = c("mock","covid")) # This will be automatically done by DESeq but it's better to specify your own levels
cell_type = c(rep("NHBE", 6),rep("A549", 6),rep("A549-ACE2", 6),rep("Calu3", 6),rep("Lung", 4))
coldata = data.frame(condition, cell_type, row.names = coldata)
head(coldata)


write.csv(coldata, file = "design_data.csv")

# Viewing distribution of reads per sample
totals = colSums(GSE_subset) #### sum counts for each column
colors <- c(rep(c("red", "hotpink", "green4", "green", "blue", "lightblue", "purple", "plum"), each = 3), rep(c("gold4", "gold"), each = 2))
barplot(totals, col = colors, names.arg = "")
legend(x = "topright", legend=c("S1_Mock","S1_Cov","S2_Mock","S2_Cov","S6_Mock","S6_Cov","S7_Mock","S7_Cov","S15_Mock","S15_Cov"), fill = c("red", "hotpink", "green4", "green", "blue", "lightblue", "purple", "plum","gold4", "gold"), bty = "n")