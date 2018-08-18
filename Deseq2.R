source("http://bioconductor.org/biocLite.R") # source of BiocLite
biocLite("DESeq2") # installing the DESeq2 package if not installed before
library(DESeq2) #load the library DESeq2
library("vsn")
biocLite('pcaExplorer') ##modules for PCA explorer
library("ggplot2") 
library("RColorBrewer")
library("gplots")
library(pcaExplorer)
library("ggbeeswarm")
getwd()
setwd("/home/tarakaramji/isomiRs/Tools/Chimira/CLLsamples/miRcounts/")## set the directory to work
counts <- read.delim("count_table.txt") #Raw counts of the samples 
#samples <- data.frame(timepoint = rep(c("ctrl", "T1", "T2"), each=5), patient=rep(1:5, 3))# for 5 patients two time points
samples <-data.frame(row.names = c('PAT1.ctrl', 'PAT2.ctrl', 'PAT3.ctrl', 'PAT4.ctrl', 'PAT5.ctrl', 'PAT1.T1', 'PAT2.T1', 'PAT3.T1','PAT4.T1', 'PAT5.T1', 'PAT1.T2', 'PAT2.T2', 'PAT3.T2','PAT4.T2', 'PAT5.T2'),timepoint = rep(c("ctrl", "T1", "T2"),each=5),patient = rep(c("PAT1","PAT2","PAT3","PAT4","PAT5")))
#ds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design=~timepoint) # matrix for DESeq and design time point

ds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design=~patient + timepoint) #matrix and design to test for differences between time points

#ds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design=~timepoint + patient) #matrix and design to test for differences between patients

ds$patient <- factor(ds$patient) # Patients (which you have as 1,2,3,4,5) is encoded as a factor and not as a numeric before running DESeq().
colnames(ds) <- colnames(counts) # include the sample names in the data set object
ds <- DESeq(ds) # DESeq2 analysis 

#Compute a scaling factor for each sample to account for differences in read depth and complexity between samples
#Estimate the variance among samples
#Test for differences in expression among groups

res_T2_T1<-results(ds, contrast=c("timepoint","T2","T1"), independentFiltering=FALSE,cooksCutoff = FALSE) # comparison of Two time points T1 and T2 into a table 
res_T2_ctrl<-results(ds, contrast=c("timepoint","T2","ctrl"), independentFiltering=FALSE,cooksCutoff = FALSE) 
res_T1_ctrl<-results(ds, contrast=c("timepoint","T1","ctrl"), independentFiltering=FALSE,cooksCutoff = FALSE) 
write.table(res_T1_ctrl,file="res_T1_ctrl", sep="\t")
write.table(res_T2_ctrl,file="res_T2_ctrl", sep="\t")
write.table(res_T2_T1,file="res_T1_T2", sep="\t")

###T1 and control comparisons
nrow(res_T1_ctrl) # display the number of rows in the result table 
sum( is.na(res_T1_ctrl$pvalue) ) # display 
res_T1_ctrl <- res_T1_ctrl[ ! is.na(res_T1_ctrl$pvalue), ]  # filtering out of Pval NA
nrow(res_T1_ctrl) # display the rows of res
sig_T1_ctrl <- res_T1_ctrl[ which(res_T1_ctrl$padj < 0.01), ] # Filtering based on Padj <0.01
sig_T1_ctrl <- sig_T1_ctrl[ order(sig_T1_ctrl$padj), ] #sorting the results by statistical significance
sig_T1_ctrl <- as.data.frame(sig_T1_ctrl) # 	table to dataframe for easy manipulation
### Options to change the display of table incase of several lines
options(width=120)  ## Display width (number of characters)
options(digits=5)   ## Number of digits to show for numbers
head(sig_T1_ctrl, n=20) # diplay first 20 lines

#####T2 and ctrl comparisons

nrow(res_T2_ctrl) # display the number of rows in the result table 
sum( is.na(res_T2_ctrl$pvalue) ) # display 
res_T2_ctrl <- res_T2_ctrl[ ! is.na(res_T2_ctrl$pvalue), ]  # filtering out of Pval NA
nrow(res_T2_ctrl) # display the rows of res
sig_T2_ctrl <- res_T2_ctrl[ which(res_T2_ctrl$padj < 0.01), ] # Filtering based on Padj <0.01
sig_T2_ctrl <- sig_T2_ctrl[ order(sig_T2_ctrl$padj), ] #sorting the results by statistical significance
sig_T2_ctrl <- as.data.frame(sig_T2_ctrl) # 	table to dataframe for easy manipulation
### Options to change the display of table incase of several lines
options(width=120)  ## Display width (number of characters)
options(digits=5)   ## Number of digits to show for numbers
head(sig_T2_ctrl, n=20) # diplay first 20 lines

#####T2 and T1 comparisons

nrow(res_T2_T1) # display the number of rows in the result table 
sum( is.na(res_T2_T1$pvalue) ) # display 
res_T2_T1 <- res_T2_T1[ ! is.na(res_T2_ctrl$pvalue), ]  # filtering out of Pval NA
nrow(res_T2_T1) # display the rows of res
sig_T2_T1 <- res_T2_T1[ which(res_T2_T1$padj < 0.01), ] # Filtering based on Padj <0.01
sig_T2_T1 <- sig_T2_T1[ order(sig_T2_T1$padj), ] #sorting the results by statistical significance
sig_T2_T1 <- as.data.frame(sig_T2_T1) # 	table to dataframe for easy manipulation
### Options to change the display of table incase of several lines
options(width=120)  ## Display width (number of characters)
options(digits=5)   ## Number of digits to show for numbers
head(sig_T2_T1, n=20) # diplay first 20 lines


###Incase you are working with Ensembl IDs. Collapse and change them to genesymbols using ENSEMBL
#source("http://bioconductor.org/biocLite.R") # biocLite souce
#biocLite("biomaRt")# if required to install
#library(biomaRt) # loading biomaRt library
#ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
#genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
 #                filters = "ensembl_gene_id",
 #                values = rownames(sig),
 #                mart = ensembl )
###data frame genemap now contains a mapping of Ensembl gene IDs to gene symbols
#head(genemap) # diplay the genemap dataframe
## Add the gene symbols to the result table
##The tapply() function call above is needed to deal with cases where there are multiple symbols for the same gene. This call maps each Ensembl gene ID to a string of one more more gene symbols separated by semi-colon
#symbols <- tapply(genemap$hgnc_symbol, genemap$ensembl_gene_id, paste, collapse="; ")
#sig$symbol <- symbols[ rownames(sig) ]
#head(sig)


###Visualizing the gene expression data
##plot visualises the differences between measurements taken in two samples, by transforming the data onto M (log ratio) and A (mean average) scales
plotMA(ds) #MA plot

## ####Apply regularized-log transform to counts
##
rld <- rlog(ds, blind=FALSE) # rlog
vsd <- varianceStabilizingTransformation(ds, blind=FALSE) # Variance stabilizing transformation very similar to rlog used with large number of samples
#vsd.fast <- vst(ds, blind=FALSE) ## fast version of Variance stabilizing transformation
#log2tr<-log2(counts(ds,normalized=TRUE)+1)# Log2 transformation with y = log2(n + 1) 
head(assay(rld), 3)  # to display first 3 rlog tranformed rows
#### Effect of transformation in the variance

notAllZero <- (rowSums(counts(ds))>0)
meanSdPlot(log2(counts(ds,normalized=TRUE)[notAllZero,] + 1)) ## log2(n+1) counts 
meanSdPlot(assay(rld[notAllZero,])) ##rlog transformation
meanSdPlot(assay(vsd[notAllZero,])) ##variance stabilizing transformation

##Smaple distance based on rlog transformation

sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$patient, rld$timepoint, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
clustering_distance_rows=sampleDists,
clustering_distance_cols=sampleDists,
col=colors)
###Plot distribution of expression values
###rld A DESeqTransform object.
###plot_type Character, choose one of boxplot, violin or density. Defaults to density
##distro_expr(rld, plot_type = "density")

rlt <- DESeq2::rlogTransformation(ds)
distro_expr(rlt,plot_type = "boxplot") ##A plot with the distribution of the expression values boxplot

## Principal component analysis
plotPCA(rld, intgroup="timepoint") # with time point 
##customize the PCA plot
data <- plotPCA(rld, intgroup=c("patient", "timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=timepoint, shape=patient)) +
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
coord_fixed()
###

pcaplot(rld,intgroup = c("patient","timepoint"),pcX = 1, pcY = 2, title = " PCA on samples - PC1 vs PC2")
pcaplot(rld, intgroup=c("timepoint"), pcX=1, pcY=4, title="PC1 vs PC2", ellipse= TRUE)  ###PCA plot on different principal components

## Heatmap of sample distances
library("gplots")   # If this fails, run: install.packages("gplots")
library("RColorBrewer")
sampleDists <- dist(t(assay(rld))) ### sample distance
sampleDistMatrix <- as.matrix( sampleDists )
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
heatmap.2(sampleDistMatrix, trace="none", col=colours)
	###removing batch effect
	assay(rld) <- limma::removeBatchEffect(assay(rld), rld$patient)
	pcaplot(rld,intgroup = c("patient","timepoint"),pcX = 1, pcY = 2, title = " PCA on samples - PC1 vs PC2")
	pcaplot(rld,intgroup = c("patient","timepoint"),pcX = 1, pcY = 3, title = " PCA on samples - PC1 vs PC3")
	pcaplot(rld,intgroup = c("patient","timepoint"),pcX = 2, pcY = 3, title = " PCA on samples - PC2 vs PC3")
############


## Heatmap of most variable genes
#library("genefilter") # loading library
#topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), 20) ### 20 variable genes
#heatmap.2(assay(rld)[topVarGenes, ], scale="row",
#trace="none", dendrogram="column", margins=c(5, 10),
#col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

###Top gene count in all the patient samples
topGene <- rownames(res_T1_ctrl)[which.min(res_T1_ctrl$padj)]
plotCounts(ds, gene = topGene, intgroup=c("patient"))

### customized top gene distribution
geneCounts <- plotCounts(ds, gene = topGene, intgroup = c("timepoint","patient"),returnData = TRUE)
ggplot(geneCounts, aes(x = timepoint, y = count, color = patient)) + scale_y_log10() +  geom_beeswarm(cex = 3)
##volcano plot

plot(res_T1_ctrl$log2FoldChange,-log(res_T1_ctrl$padj),points(sig_T1_ctrl$log2FoldChange, -log(sig_T1_ctrl$padj),col="red"))

####Clear list and restart
rm(list = ls())
.rs.restartR()
