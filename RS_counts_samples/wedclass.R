## Set your working directory
setwd("~/Desktop/Spring 2020/SCHOOL/Genomics Class/RS_counts_samples")

## Import the libraries that we're likely to need in this session
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  ### First: BiocManager::install("vsn") AND BiocManager::install("hexbin")

## Import the counts matrix
countsTable <- read.table("RS_cds2kb_countsMatrix.txt", header=TRUE, row.names=1)

head(countsTable)
dim(countsTable)
countsTableRound <- round(countsTable) # Need to round because DESeq wants only integers
head(countsTableRound)

## Import the samples description table - links each sample to factors of the experimental design.
# Need the colClasses otherwise imports "day" as numeric which DESeq doesn't like, coula altneratively change to d0, d5, d10
conds <- read.delim("RS_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1, colClasses=c('factor', 'factor', 'factor', 'factor'))
head(conds)
dim(conds)


## Let's see how many reads we have from each sample:
colSums(countsTableRound)
mean(colSums(countsTableRound))
barplot(colSums(countsTableRound), las=3, cex.names=0.5,names.arg = substring(colnames(countsTableRound),1,13))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd =2)
#average number of counts per gene
mean(rowSums(countsTableRound))
median(rowSums(countsTableRound))
# dispersion across genes - differences 

# avergae number of counts per gene per sample
apply(countsTableRound,2,mean)

## Create a DESeq object and define the experimental design here with the tilde
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, design = ~ climate + day + treatment)

dim(dds)
# [1] 66408    76
# Filter out genes with few reads

dds <- dds[rowSums(counts(dds))>760]
dim(dds)
# [1] 23887    76 Filterring to sum of 76 reads across all samples

dds <- DESeq(dds)

#List the results you've generated

resultsNames(dds)

#Running the model with pop
# [1] "Intercept"            "pop_BRU_05_vs_ASC_06"
# [3] "pop_CAM_02_vs_ASC_06" "pop_ESC_01_vs_ASC_06"
# [5] "pop_JAY_02_vs_ASC_06" "pop_KAN_04_vs_ASC_06"
# [7] "pop_LOL_02_vs_ASC_06" "pop_MMF_13_vs_ASC_06"
# [9] "pop_NOR_02_vs_ASC_06" "pop_XBM_07_vs_ASC_06"
# [11] "day_10_vs_0"          "day_5_vs_0"          
# [13] "treatment_D_vs_C"     "treatment_H_vs_C" 

# Running the model with climate
# [1] "Intercept"        "climate_HD_vs_CW" "day_10_vs_0"  
# [4] "day_5_vs_0"       "treatment_D_vs_C" "treatment_H_vs_C"


## Run the DESeq model to test for differential gene expression: 1) estimate size factors (per sample), 2) estimate dispersion (per gene), 3) run negative binomial glm



# List the results you've generated



# Order and list and summarize results from specific contrasts
# Here you set your adjusted p-value cutoff, can make summary tables of the number of genes differentially expressed (up- or down-regulated) for each contrast
res <- results(dds, alpha = 0.05)
res <- res[order(res$padj),]
head(res)


summary(res)

res_treatCD <- results(dds, name = "treatment_D_vs_C", alpha=0.05)

res_treatCD <- results(dds, alpha = 0.05)
res_treatCD <- res[order(res$padj),]
head(res_treatCD)

# log2 fold change (MLE): treatment H vs C 
# Wald test p-value: treatment H vs C 
# DataFrame with 6 rows and 6 columns
# baseMean    log2FoldChange             lfcSE
# <numeric>         <numeric>         <numeric>
#   MA_172878g0010   15.8548874481417  2.26913549236445  0.43955520945096
# MA_28973g0010    18.8813749792546 -1.96620213947109 0.411529496517035
# MA_10426002g0010 10.8980752578363 -1.20714151942448 0.281202070702332
# MA_10426033g0010 67.6992372010061 0.738619854758194 0.182127301545519
# MA_10429525g0010 60.5937645508928  1.17170143286339 0.281637964070113
# MA_10433965g0010 49.9866944326994 0.612714461794651 0.151941670249246
# stat               pvalue
# <numeric>            <numeric>
#   MA_172878g0010    5.16234466928234 2.43875719941453e-07
# MA_28973g0010    -4.77779152190056  1.7723099296183e-06
# MA_10426002g0010 -4.29279029279376  1.7644163823468e-05
# MA_10426033g0010   4.0555141842564 5.00241365216353e-05
# MA_10429525g0010  4.16031069082609 3.17814956479126e-05
# MA_10433965g0010  4.03256368571934 5.51716566884857e-05
# padj
# <numeric>
#   MA_172878g0010   0.00147447260276602
# MA_28973g0010    0.00535769291723613
# MA_10426002g0010  0.0256590643337373
# MA_10426033g0010  0.0256590643337373
# MA_10429525g0010  0.0256590643337373
# MA_10433965g0010  0.0256590643337373

summary(res_treatCD)


##### Data visualization #####
# MA plot
plotMA(res_treatCD,ylim=c(-3,3))

# PCA
vsd <- vst(dds, blind=FALSE)

data <- plotPCA(vsd, intgroup=c("climate","treatment","day"),returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
data$treatment <- factor(data$treatment, levels=c("C","H","D"), labels = c("C","H","D"))
data$day <- factor(data$day, levels=c("0","5","10"), labels = c("0","5","10"))

ggplot(data, aes(PC1, PC2, color=climate, shape=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()


  
# Counts of specific top gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="MA_10257300g0010", intgroup = (c("treatment","climate","day")), returnData=TRUE)
d



p <-ggplot(d, aes(x=treatment, y=count, color=day, shape=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
  scale_x_discrete(limits=c("C","H","D"))
p


# Heatmap of top 20 genes sorted by pvalue
library(pheatmap)
topgenes <- head(rownames(res_treatCD),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("treatment","climate","day")])
pheatmap(mat, annotation_col=df)



