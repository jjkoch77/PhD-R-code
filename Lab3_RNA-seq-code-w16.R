####################################################################################
##
## Bioinf545 / Stats545 / Biostat646
##  Lab 3 - Testing differential expression with RNA-seq data
##
####################################################################################

####  In R:
getwd()
setwd(...)

###  Input the gene count data
gene.cts<-read.table("HPV_read_count.txt",header=T,sep="\t",stringsAsFactors=F)
# Use the line below if you are working on bcs2 class-XX machine
#gene.cts<-read.table("/class/data/bio545w16/lab3/HPV_read_count.txt", header=TRUE,
 sep="\t", stringsAsFactors=FALSE)

meta<- read.table("HPV_metadata.txt",header=T,sep="\t",stringsAsFactors=F)
# Use the line below if you are working on bcs2 class-XX machine
#meta<- read.table("/class/data/bio545w16/lab3/HPV_metadata.txt",header=T,sep="\t",stringsAsFactors=F)
meta
hpv.index<- match(colnames(gene.cts[,2:19]),meta$HPV.ID)
meta<-meta[hpv.index,]
meta

gene.cts[1:4,]
dim(gene.cts)

####  Look at properties of the gene count data
summary(gene.cts[,3])
avg.raw <- colMeans(gene.cts[,2:19])
barplot(avg.raw)
#colSums(gene.cts[,2:19])


##############################################################################
##  Using the edgeR package to test differential expression
##############################################################################
library(edgeR)

## Create the design matrix
design<- cbind(int=1, subtype2.1=as.factor(meta[,2]))
design

subtype <- factor(meta[,2])
design<- model.matrix(~1+subtype)
colnames(design)<-c('int','subtype2.1')
design

## Prepare data in right format for edgeR
counts <- gene.cts[,2:19]
rownames(counts)<-gene.cts[,1]
y<-DGEList(counts=counts,group=subtype)
y<-calcNormFactors(y)  # Normalize using TMM method
y$samples # norm factors are in last column


##  Estimate a common dispersion factor (use the same for each gene) - 
##     not recommended unless the tagwise dispersion factors are similar and there's minimal trend
y<-estimateCommonDisp(y, rowsum.filter=5)
attributes(y)
y$pseudo.counts[1:3,]  # normalized counts
y$common.dispersion

##  Assess correlations among samples
library(gplots)
corrs.y<- cor(log2(y$pseudo.counts+1))
colnames(corrs.y) <- meta[,1]
rownames(corrs.y) <- meta$subtype
heatmap.2(corrs.y, symm=T, main="sample correlations", margins=c(7,7)) #,dendrogram="none")

###  "common" vs "trended" vs "tagwise" dispersion estimation
#y <- estimateDisp(y, design)

##### Fit the negative binomial glm (common dispersion) in edgeR #####
#  We're skipping this because not recommended in this case.
#  May want to use for experiments with very small sample size (when trend may not be accurately estimated)
#fit<-glmFit(y,design, dispersion=y$common.dispersion)
# Results from the above line of code include coefficients and residual degrees of freedom, but not p-values
#lrt.common<-glmLRT(fit,coef=2) # Calculate p-values for the 2nd column of the design matrix
#attributes(lrt.common)
#lrt.common$dispersion
#lrt.common$table[1:4,]
#pvals_c <- lrt.common$table$PValue
#bmp("pvals-commonDisp.bmp",width=400,height=400)
#hist(pvals_c)
#dev.off()
#bmp("pvals-commonDisp-1.bmp",width=400,height=400)
#hist(pvals_c[rowMeans(counts)>=1]) # Spike at p-value near 1 is due to 
#genes with extremely low counts
#dev.off()
#topTags(lrt.common)


### Fit model with trended dispersion (as a function of expression level or abundance)
y<-estimateGLMTrendedDisp(y,design)
fit_tr<-glmFit(y,design, dispersion=y$trended.dispersion)
lrt.trend<-glmLRT(fit_tr,coef=2) # coefficient=2, b/c we want the p-value for the 2nd column of the design matrix.

pvals_tr <- lrt.trend$table$PValue
hist(pvals_tr)
#bmp("pvals-trended-1.bmp",width=400,height=400)
hist(pvals_tr[rowMeans(counts)>=1])
#dev.off()
topTags(lrt.trend)


###  Fit model with tagwise dispersion (separately for each gene, and then uses
##  empirical Bayesian model to "shrink" estimates toward trended (or common) overdispersion estimate.
y<-estimateGLMTagwiseDisp(y,design,trend=TRUE) #If trend=False, then the trended estimates are not used for the background
fit_tag<-glmFit(y,design)
lrt.tagwise<-glmLRT(fit_tag,coef=2)

pvals_tag <- lrt.tagwise$table$PValue
hist(pvals_tag)
hist(pvals_tag[rowMeans(counts)>=1])
topTags(lrt.tagwise)

# Visualize the relationship between average logCPM and dispersion estimate
#bmp("AveLogCPM_vs_dispersion.bmp",width=300,height=300)
plot(y$AveLogCPM,y$tagwise.dispersion, ylim=c(0,6),pch=".")
lines(sort(y$AveLogCPM),y$trended.dispersion[order(y$AveLogCPM)], col="red",lwd=2) # trend line
#dev.off()

## Volcano plot
logFC_tag <- lrt.tagwise$table$logFC
#bmp("VolcanoPlot-edgeR.bmp",width=300,height=300)
plot(logFC_tag,-log10(pvals_tag),xlab="Log2 fold change",ylab="-log10(p-values)",
main="Volcano plot of tagwise results",xlim=c(-15,15))
#dev.off()

##  Calculate False Discovery Rate (FDR) 
FDR_tag<- p.adjust(pvals_tag, method="BH") # "BH" stands for the Benjamini-Hochberg FDR method
FDRtag_no0s<-p.adjust(pvals_tag[rowMeans(counts)>=1])  # may want to filter out gene with very low count first


## Merge in Entrez gene IDs 
library(org.Hs.eg.db)
xx <- as.list(org.Hs.egALIAS2EG)
# Remove identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]
## Get the first one, when there are more than 1 possible entrez ID
entrez.all<-unlist(lapply(xx,function(yy){ as.numeric(yy[[1]])}))
entrez.index<-match(gene.cts$symbol,names(entrez.all))
entrezid<-entrez.all[entrez.index]
entrezid[1:10]
length(entrezid[is.na(entrezid)])  # Shows how many genes did not have any matched Entrez ID.

## Output tagwise results to a file
to.output<-cbind(entrezid,gene.cts,lrt.tagwise$table,FDR_tag)
to.output[1:4,]
write.table(to.output, file="Lab3-RNAseq_DE_results.txt", sep="\t", col.names=T, row.names=F)

##  Comparison of trended versus tagwise edgeR results
#bmp("edgeR-pvalue-comparison.bmp",width=600,height=600)
par(mfrow=c(1,2))
plot(pvals_tr, pvals_tag, pch=".")
plot(-log10(pvals_tr), -log10(pvals_tag), pch=".")

#plot(pvals_c, pvals_tag, pch=".")
#plot(-log10(pvals_c), -log10(pvals_tag), pch=".")
#dev.off()


##############################################################################
##   Use DESeq2 R package to test for differentially expressed genes
##############################################################################

library(DESeq2)
#clean.cts<-as.matrix(round(gene.cts[,2:19])) # Use IF not all data were integers
coldata<-data.frame(group=subtype)
cds<-DESeqDataSetFromMatrix(gene.cts[,2:19],colData=coldata,design=~group)
#cds2<-DESeqDataSetFromMatrix(gene.cts[,2:19],colData=coldata,design=~1+group)
# Look at size factors (can skip this, because it is done in deseq function)
rownames(cds)<-gene.cts[,1] 
cds <- estimateSizeFactors(cds)
sizeFactors(cds) 

##########  Run DEseq2 Analysis  ###########
cds <- DESeq(cds) # Does the Empirical Bayes test
deseq2.res <- results(cds)  # Extract results table
deseq2.res<-deseq2.res[order(deseq2.res$pvalue),]
deseq2.res[1:6,]
summary(deseq2.res,alpha=0.05)

###  DESeq2  MA-plot  ###
#bmp("DESeq2-MAplot.bmp",width=400,height=400)
plotMA(deseq2.res, main="DESeq2", ylim=c(-10,8))
#dev.off()

###  DESeq2 - Expression versus standard error  ###
#bmp("meanExpr_vs_SE-DESeq2.bmp",width=400,height=400)
plot(log2(deseq2.res$baseMean),deseq2.res$lfcSE,xlab="baseMean expression",
 ylab="Standard Error",ylim=c(0,1.1),cex=0.25)
grid()
#dev.off()

#########  Histogram of DESeq2 p-values  ##########
#bmp("pvalue-histogram-DESeq2.bmp",width=300,height=300)
hist(deseq2.res$pvalue)
#dev.off()


## Compare -log10 pvalues
unsort.res<-results(cds) # Get unsorted DESeq2 results
logEdger<-(-1)*log10(pvals_tag+10^-300) # a small constant is added to avoid taking log of zero
logdeseq<-(-1)*log10(unsort.res$pval+10^-300)

#bmp("Compare_qqPlots_RNAseq.bmp",width=650,height=250)
par(mfrow=c(1,3))
plot(logEdger,logdeseq,main="All genes",cex=0.25)
lines(c(0,150),c(0,150),col="red")
plot(logEdger[rowMeans(counts)>=500],logdeseq[rowMeans(counts)>=500],main="High genes",cex=0.25)
lines(c(0,150),c(0,150),col="red")
plot(logEdger[rowMeans(counts)<500],logdeseq[rowMeans(counts)<500],main="Low genes",cex=0.25)
lines(c(0,150),c(0,150),col="red")
#dev.off()

##############################################################################
##   Use limma voom to test for differentially expressed genes
##############################################################################
design<- model.matrix(~1+subtype)
colnames(design)<-c('int','subtype2.1')
design

####  The voom transformation  ####
##  This models the mean-variance relationship and calculates per-observation weights that will be carried over to the lmFit and eBayes functions of limma
voom.y <- voom(counts=y, design=design,plot=T)
voom.yq <- voom(y, design, normalize.method= "quantile", plot=T) 

####  Perform linear analysis (t-test)  ####
fit <- lmFit(voom.y, design)
attributes(fit)  # Look at what’s in the results
fit$coefficients[1:4,]  # 1st column is intercept, 2nd column are log2 fold changes.
fit$sigma[1:4]  # sample standard deviations
fit$df.residual[1:4]  # degrees of freedom for a standard t-test
fit$stdev.unscaled[1:4,] # multiplier for the st dev to get denominator for t-statistic

denom <- fit$stdev.unscaled[,2]*fit$sigma
tstat <- fit$coefficients[,2]/denom
pvalue <- 2*pt(abs(tstat),fit$df.residual,lower.tail=F)

####   Use the empirical Bayes Test   ####
# Test for differential expression (empirical Bayes - eBayes)
fit <- eBayes(fit)  
attributes(fit)  # see what’s in the results
head(fit$p.value)  # look at top significant p-values
fdr.eBayes <- p.adjust(fit$p.value[,2], method="BH") # adjust p-values for multiple comparisons
hist(fit$p.value[,2])
sort(fdr.eBayes)[1:4]

####   Calculate fold changes   ####
log2fold<-fit$coefficients[,2]
# Two possible ways to calculate fold changes
foldchg1 <- 2^log2fold
range(foldchg1)  # In this one, down-regulated genes are [0,1] while up-regulated genes are in (1,infinity)
foldchg2 <- ifelse(log2fold>0, 2^log2fold,-1/(2^log2fold))
range(foldchg2)  # In this one, up-regulated genes are the same as above, while down-regulated genes are flipped over to the negative range.
# For example, if a gene is down-regulated from a value of 100 down to 25, 
# the first method would report its  fold change as 0.25, while the second would report -4 fold.



##############################################################################################
##  For illustrative purposes only, this shows how to use Poisson regression to test differential expression
##   Not recommended, because usually this greatly underestimates the variance

##  First, just test one gene to see how it works
clean.cts[100,]
fit1<- glm(as.numeric(clean.cts[100,]) ~ factor(group), family = poisson)
G1<-summary(fit1)
attributes(G1)  # We want coefficients
G1$coefficients
G1$coefficients[2,1]  # estimate
G1$coefficients[2,4]  # p-value of interest

##  Now test all genes in a loop
length(clean.cts[,1])  # 33159
p.est<-NA
p.pval<-NA
for (i in 1:33159) {
	temp<-glm(as.numeric(clean.cts[i,])~factor(group), family=poisson)
	g1<-summary(temp)
	p.est[i]<-g1$coefficients[2,1]
	p.pval[i]<-g1$coefficients[2,4]
}
# The above loop may take a couple minutes to run

