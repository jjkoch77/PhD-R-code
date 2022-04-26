##Running csaw to examine the 5-hydroxy methylation results (HMeDIP-seq)

##Install dependency (csaw and edgeR packages):
source("https://bioconductor.org/biocLite.R")
biocLite("csaw")
biocLite("edgeR")
#First, load in R packages:
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(csaw)
library(edgeR)
library(rtracklayer)

############################################################################################

##Note: BAM files for HMeDIP-seq data are quite large (~4 GB per sample), so run this on the ONES server!

##Set working directory
getwd()
#ONES Server working directory: 
setwd('/jjkoch/projects/dolinoydriftseq/pulldown/bowtie2_bams/')

##################################################################################################
##1. Loading in data from BAM files
#Example from real data:
# bam.files <- c("es_1.bam", "es_2.bam", "tn_1.bam", "tn_2.bam")
# design <- model.matrix(~factor(c("es", "es", "tn", "tn")))
# colnames(design) <- c("intercept", "cell.type")

#Set up multicore processing to speed up R code:
#multicoreParam <- MulticoreParam(workers = 4)

#restrict.param <- readParam(restrict=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
#"chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
#"chr18","chr19"), BPPARAM=multicoreParam
# remove sex chromosomes from analysis and use multiple cores for processing
restrict.param <- readParam(restrict=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                                       "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
                                       "chr18","chr19"))

#ALL 2-4-10M samples:
bam.files = c('101a_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam','241c_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '431a_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam','441d_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '951a_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam','991d_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '101a_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam','241c_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '431a_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam','441d_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '951a_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam','991d_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '101a_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam','241c_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '431a_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam','441d_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '951a_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam','991d_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '161a_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam','171c_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '501a_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam','411b_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '941a_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam','981c_2M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '161a_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam','171c_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '501a_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam','411b_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '941a_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam','981c_4M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '161a_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam','171c_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '501a_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam','411b_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam',
              '941a_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam','981c_10M_hmc_pulldown_trimmed.fq.gz_aligned.bam')

data <- windowCounts(bam.files, ext=52, width=100, param=restrict.param)
#Note: Tried to do all 24 samples combined on the ONES server, but the connection kept resetting prior to the function finishing;
##Turned on the SSH keep alive setting in MobaXterm (under settings --> SSH) (3/27/17)

#Lab computer error on 3/9/17: [W::sam_hdr_read] bgzf_check_EOF: Invalid argument
##3/10/17 -- No error on the ONES server!

#Note: window refers to the expected peak size. For TF binding, this will be small,
#but for histone marks, it will be larger (>/= 150 bp).
#Given that 5hmC is a diffuse mark, the window size should err on the side of small;
#here, I've put down 20 bp for the window width.

#ext refers to fragment length of directional read extension.

#Testing out smaller window width (1/8/18):
data <- windowCounts(bam.files, ext=52, width=20, param=restrict.param)



##################################################################################################
##2. Filter out uninteresting regions

#Note: It is difficult to know where to draw the line. Start by choosing a strategy:
#I have chosen Strategy #4 in the csaw vignette - filter by local enrichment (tiled background estimates for each window)
#This method is similar to MACS, which is used in the mint pipeline:
surrounds <- 2000
neighbor <- suppressWarnings(resize(rowRanges(data), surrounds, fix="center"))
wider <- regionCounts(bam.files, regions=neighbor, ext=52, param=restrict.param)
# Warning message in "wider" step (see below) on 3/28/17:Warning message:
#In valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE) :
#  GRanges object contains 40 out-of-bound ranges located on sequence
#chrM. Note that only ranges located on a non-circular sequence whose
#length is not NA can be considered out-of-bound (use seqlengths() and
#isCircular() to get the lengths and circularity flags of the underlying
#sequences). You can use trim() to trim these ranges. See
#?`trim,GenomicRanges-method` for more information.
filter.stat <- filterWindows(data, wider, type="local")
#keep1 <- filter.stat$filter > log2(3) #4/4/17 -- 3 fold enrichment over surrounding background required for inclusion
keep1 <- filter.stat$filter > log2(2) #4/10/17 -- 2 fold enrichment over surrounding background required for inclusion
sum(keep1) #3/31/17 -- all control samples: [1] 1626190
#4/4/17 -- all 2M-10M samples: [1] 1522880
#4/10/17 -- all 2M-10M samples (2 fold change): [1] 3635436
#5/18/17 - full dataset -- [1] 3629322

##Count size filter to further remove spurious results:
keep2 <- keep1 > aveLogCPM(5, lib.size=mean(data$totals))
sum(keep2) #3/31/17 -- all control samples: [1] 47132875
#4/4/17 -- all 2M-10M samples: [1] 52457553
#5/18/17 - [1] 47898974
#Note: Different filters can also be combined in more advanced applications, e.g., by running
filtered.data <- data[keep1 & keep2,] #for filter vectors keep1 and keep2
##Combine local background and count-based filters for more accurate data! Want to ensure
#local enrichment over background still has a high count!

############################################################################################
##3. Calculate normalization factors
##Note: PePr is using a modified TMM method to normalize for the difference in IP efficiencies between samples 
#It is making an implicit assumption that there is substantial overlap of peaks in every sample. 
#However, it is sometimes not true between groups (for example, between TF ChIP-seq and TF knockout). 
#So for differential binding analysis, switch to intra-group normalization.
##In csaw, binned (composition) normalization corresponds to the intra-group normalization
binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=restrict.param)
normfacs <- normOffsets(binned)
#normfacs 
#[1] 0.9774151 0.9739487 0.9912091 0.9925577 1.0127311 1.0084030 0.9950691
#[8] 0.9906483 1.0125539 1.0041050 0.9440516 1.0395918 0.9792655 1.0008891
#[15] 1.0070468 1.0248470 0.9892408 0.9578531 0.9757565 1.0242452 0.9858873
#[22] 0.9649796 0.9756632 1.0737917 0.9648474 0.9695977 1.0193505 0.9828881
#[29] 0.9688325 1.0323669 0.9841850 1.0236111 1.1749561 1.1618832 0.9577115
#[36] 0.9015898

#Note: Generally, they are all quite close already, which is good!

##Note: bin=10000 is appropriate for most datasets, but test multiple bin sizes!
#Examples from vignette:
#demo <- windowCounts(bam.files, bin=TRUE, width=5000, param=param)
#normOffsets(demo)
#[1] 1.0093596 0.9754046 1.0105468 1.0051082
#demo <- windowCounts(bam.files, bin=TRUE, width=15000, param=param)
#normOffsets(demo)
#[1] 1.007559 0.972134 1.016096 1.004775
##When these values are different from each other, normalization is required!

##3/31/17 -- MA plot line goes straight through the center of the clouds. Good!

############################################################################################
##4. Identifying regions of differential peaks (DHMRs)
#Establish experimental design
#From edgeR package vignette -- The full interaction formula is
#design <- model.matrix(~Treat + Time + Treat:Time, data=targets)
#fit <- glmFit(y, design)
#This formula is primarily useful as a way to conduct an overall test for interaction. 
#The coefficient names are:
#colnames(design)
#[1] "(Intercept)" "TreatDrug"
#[3] "Time1h" "Time2h"
#[5] "TreatDrug:Time1h" "TreatDrug:Time2h"

#For the full dataset, our design is as follows:
n <- c(1:36)
age <- c("2M", "2M", "2M", "2M", "2M", "2M","4M", "4M", "4M", "4M", "4M", "4M","10M", "10M", "10M", "10M", "10M", "10M",
         "2M", "2M", "2M", "2M", "2M", "2M","4M", "4M", "4M", "4M", "4M", "4M","10M", "10M", "10M", "10M", "10M", "10M")
exposure <- c("Control", "Control","Control", "Control","Control", "Control","Control", "Control","Control", "Control","Control", "Control","Control", "Control","Control", "Control","Control", "Control",
              "BPA","BPA","BPA","BPA","BPA","BPA","BPA","BPA","BPA","BPA","BPA","BPA","BPA","BPA","BPA","BPA","BPA","BPA")
mouse <-c(1:6,1:6,1:6,7:12,7:12,7:12) #paired mouse covariate
sex <- c("F","M", "F","M", "F","M","F","M", "F","M", "F","M","F","M", "F","M", "F","M", 
         "F","M", "F","M", "F","M","F","M", "F","M", "F","M","F","M", "F","M", "F","M")
design_combined <- data.frame(n,age,exposure,mouse,sex)
design <- model.matrix(~age + exposure + age:exposure + mouse + sex, data=design_combined)

##Setting up the data:
y <- asDGEList(filtered.data, norm.factors=normfacs)

#Empirical bayes method used to stabilize estimates
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.01399 0.01551 0.01797 0.02048 0.02283 0.05584

#4/4/17: Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     0.04074 0.04217 0.04559 0.04684 0.05039 0.06942

#4/10/17    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#       0.03172 0.03426 0.04062 0.04325 0.05088 0.06857


##Fit model for differential 5hmC binding; Empirical bayesian method used to 
#stabilize the QL dispersion estimates
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.post)
head(fit$coefficients)

results_age <- glmQLFTest(fit, coef="age2M") #DB by age (10M to 2M; 2M is reference!!) is being tested here
results_exposure <- glmQLFTest(fit, coef="exposureControl") #DB by exposure (BPA to Control; Control is reference!!) is being tested here
results_sex <- glmQLFTest(fit, coef="sexM") #DB by sex (Female to Male; M is reference!!) is being tested here
results_age_exp <- glmQLFTest(fit, coef="age2M:exposureControl") #DB by age:exp interaction is being tested here

#The null hypothesis here is that the cell type has no effect. The contrast argument in
#the glmQLFTest function specifies which factors are of interest. In this case, a contrast of
#c(0, 1) defines the null hypothesis as 0*intercept + 1*cell.type = 0, i.e., that the log-fold
#change between cell types is zero. DB windows can then be identified by rejecting the null.

#Note: can also visualize effects of empirical bayes method on dispersion (SEE VIGNETTE).

############################################################################################
##5. Correct for multiple testing
#Cluster the samples (quick and dirty method) using defined window length (tol)
#Note:The chosen tol represents the minimum distance at which two binding events are treated as separate sites.
#Large values (500 - 1000 bp) reduce redundancy and favor a region-based interpretation of
#the results, while smaller values (< 200 bp) allow resolution of individual binding sites.
merged <- mergeWindows(rowRanges(filtered.data), tol=500L) #932609 ranges
merged$region

#merged <- mergeWindows(rowRanges(filtered.data), tol=200L) #551229  ranges
#merged$region

#A simple check can be used to determine whether most clusters are of an
#acceptable size. Huge clusters indicate that more aggressive filtering from Chapter 3 is
#required. This mitigates chaining effects by reducing the density of windows in the genome.
summary(width(merged$region))
#4/4/17:    
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#100.0   150.0   200.0   203.4   250.0  1500.0

#4/10/17 AND 5/22/17:
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#100     200     250     274     300    3150

#Combined p-value is computed for each cluster with the combineTests function; The
#BH method is then applied to control the FDR across all detected clusters.
tabcom_age <- combineTests(merged$id, results_age$table)
#9/15/17 - write.csv(tabcom_age, file = "DHMRs_csaw_tabcom_age_091517.csv")

tabcom_exposure <- combineTests(merged$id, results_exposure$table)
#9/15/17 - write.csv(tabcom_exposure, file = "DHMRs_csaw_tabcom_exposure_091517.csv")

tabcom_sex <- combineTests(merged$id, results_sex$table)
#9/15/17 - write.csv(tabcom_sex, file = "DHMRs_csaw_tabcom_sex_091517.csv")

tabcom_age_exp <- combineTests(merged$id, results_age_exp$table)
#9/15/17 - write.csv(tabcom_age_exp, file = "DHMRs_csaw_tabcom_age.exposure_091517.csv")


head(tabcom_age)
head(tabcom_exposure)
head(tabcom_sex)
head(tabcom_age_exp)
is.sig.age <- tabcom_age$FDR <= 0.10
is.sig.exp <- tabcom_exposure$FDR <= 0.10
is.sig.sex <- tabcom_sex$FDR <= 0.10
is.sig.age.exp <- tabcom_age_exp$FDR <= 0.10
has.up.age <- tabcom_age$logFC.up > 0
has.up.exp <- tabcom_exposure$logFC.up > 0
has.up.sex <- tabcom_sex$logFC.up > 0
has.up.age.exp <- tabcom_age_exp$logFC.up > 0
has.down.age <- tabcom_age$logFC.down > 0
has.down.exp <- tabcom_exposure$logFC.down > 0
has.down.sex <- tabcom_sex$logFC.down > 0
has.down.age.exp <- tabcom_age_exp$logFC.down > 0

#Breakdown by model coefficients (4/4/17):
dhmrs_age <- data.frame(Total.DB=sum(is.sig.age), DB.up=sum(is.sig.age & has.up.age & !has.down.age),
           DB.down=sum(is.sig.age & !has.up.age & has.down.age), DB.both=sum(is.sig.age & has.up.age & has.down.age))
dhmrs_age

#4/4/17:
#  Total.DB DB.up DB.down DB.both
#1      239   236       3       0

#4/10/17:
#  Total.DB DB.up DB.down DB.both
#1     5368  5292      72       4

#5/22/17
#  Total.DB DB.up DB.down DB.both
#1     8613  8446     160       7

dhmrs_exposure <- data.frame(Total.DB=sum(is.sig.exp), DB.up=sum(is.sig.exp & has.up.exp & !has.down.exp),
           DB.down=sum(is.sig.exp & !has.up.exp & has.down.exp), DB.both=sum(is.sig.exp & has.up.exp & has.down.exp))
dhmrs_exposure
#4/4/17:
#  Total.DB DB.up DB.down DB.both
#1      670   328     340       2

#4/10/17:
#  Total.DB DB.up DB.down DB.both
#1     1284   624     612      48

#5/22/17
#  Total.DB DB.up DB.down DB.both
#1     5950  4247    1559     144


dhmrs_sex <- data.frame(Total.DB=sum(is.sig.sex), DB.up=sum(is.sig.sex & has.up.sex & !has.down.sex),
           DB.down=sum(is.sig.sex & !has.up.sex & has.down.sex), DB.both=sum(is.sig.sex & has.up.sex & has.down.sex))
dhmrs_sex
#4/4/17 (chrX and chrY included)
#  Total.DB DB.up DB.down DB.both
#1    33482 11767   21671      17

#4/10/17 (no chrX and chrY):
#  Total.DB DB.up DB.down DB.both
#1     2034  2007      24       3

#5/22/17 (No ChrX and chrY):
#Total.DB DB.up DB.down DB.both
#1     2944  2834     105       2


dhmrs_age.exp <- data.frame(Total.DB=sum(is.sig.age.exp), DB.up=sum(is.sig.age.exp & has.up.age.exp & !has.down.age.exp),
           DB.down=sum(is.sig.age.exp & !has.up.age.exp & has.down.age.exp), DB.both=sum(is.sig.age.exp & has.up.age.exp & has.down.age.exp))
dhmrs_age.exp
#4/4/17 AND 4/10/17 AND 5/22/17:
#  Total.DB DB.up DB.down DB.both 
#1        0     0       0       0

############################################################################################
##6. Annotate the called regions!

#Next, create annotatios of merged regions from the UCSC mm10 reference genome
anno <- detailRanges(merged$region, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(3000, 1000), dist=5000) #promoter defined as 3000 bp upstream to 1000 bp downstream
#Examine top sites
head(anno$overlap) #Each pattern contains GENE|EXONS|STRAND to describe the strand and overlapped exons of that gene
#Promoters are labelled as exon 0 whereas introns are labelled as I.
head(anno$left) #For left and right, an additional [DISTANCE] field is included. 
#This indicates the gap between the annotated feature and the supplied region.
head(anno$right)

#Note: While the string representation saves space in the output, it is not easy to work with.
#If the annotation needs to manipulated directly, users can obtain it from the detailRanges
#command by not specifying the regions of interest. This can then be used for interactive
#manipulation, e.g., to identify all genes where the promoter contains DB sites.
anno.ranges <- detailRanges(txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, orgdb=org.Mm.eg.db)
anno.ranges #GRanges object for all regions

##Save the file in a BED format using rtracklayer:

#5/22/17 -- ALL 36 samples (2,4,10 months):
is.sig.age <- tabcom_age$FDR <= 0.10
test.age <- merged$region[is.sig.age]
test.age$score <- -10*log10(tabcom_age$FDR[is.sig.age])
names(test.age) <- paste0("region", 1:sum(is.sig.age))
export(test.age, "Differential 5-hmC peaks_2-4-10M_Age_ALLSamples_FDR<0.1_052217.bed")
head(read.table("Differential 5-hmC peaks_2-4-10M_Age_ALLSamples_FDR<0.1_052217.bed"))

is.sig.exp <- tabcom_exposure$FDR <= 0.10
test.exp <- merged$region[is.sig.exp]
test.exp$score <- -10*log10(tabcom_exposure$FDR[is.sig.exp])
names(test.exp) <- paste0("region", 1:sum(is.sig.exp))
export(test.exp, "Differential 5-hmC peaks_2-4-10M_Exposure_ALLSamples_FDR<0.1_052217.bed")
head(read.table("Differential 5-hmC peaks_2-4-10M_Exposure_ALLSamples_FDR<0.1_052217.bed"))

is.sig.sex <- tabcom_sex$FDR <= 0.10
test.sex <- merged$region[is.sig.sex]
test.sex$score <- -10*log10(tabcom_sex$FDR[is.sig.sex])
names(test.sex) <- paste0("region", 1:sum(is.sig.sex))
export(test.sex, "Differential 5-hmC peaks_2-4-10M_Sex_ALLSamples_FDR<0.1_052217.bed")
head(read.table("Differential 5-hmC peaks_2-4-10M_Sex_ALLSamples_FDR<0.1_052217.bed"))

is.sig.age.exp <- tabcom_age_exp$FDR <= 0.10 # No results; doesn't work! (4/10/17 AND 5/22/17)
test.age.exp <- merged$region[is.sig.age.exp]
test.age.exp$score <- -10*log10(tabcom_age_exp$FDR[is.sig.age.exp])
names(test.age.exp) <- paste0("region", 1:sum(is.sig.age.exp))
export(test.age.exp, "Differential 5-hmC peaks_2-4-10M_Age:Exp_ALLSamples_FDR<0.1_052217.bed")
head(read.table("Differential 5-hmC peaks_2-4-10M_Age:Exp_ALLSamples_FDR<0.1_052217.bed"))

#9/19/17 -- Hypermethylated sites by age for ALL 36 samples (2,4,10 months):
is.sig.age <- tabcom_age$FDR <= 0.10
test.age <- merged$region[is.sig.age]
test.age$score <- -10*log10(tabcom_age$FDR[is.sig.age])
names(test.age) <- paste0("region", 1:sum(is.sig.age))
export(test.age, "Differential 5-hmC peaks_2-4-10M_Age_ALLSamples_FDR<0.1_052217.bed")
head(read.table("Differential 5-hmC peaks_2-4-10M_Age_ALLSamples_FDR<0.1_052217.bed"))

is.sig.exp <- tabcom_exposure$FDR <= 0.10
test.exp <- merged$region[is.sig.exp]
test.exp$score <- -10*log10(tabcom_exposure$FDR[is.sig.exp])
names(test.exp) <- paste0("region", 1:sum(is.sig.exp))
export(test.exp, "Differential 5-hmC peaks_2-4-10M_Exposure_ALLSamples_FDR<0.1_052217.bed")
head(read.table("Differential 5-hmC peaks_2-4-10M_Exposure_ALLSamples_FDR<0.1_052217.bed"))

is.sig.sex <- tabcom_sex$FDR <= 0.10
test.sex <- merged$region[is.sig.sex]
test.sex$score <- -10*log10(tabcom_sex$FDR[is.sig.sex])
names(test.sex) <- paste0("region", 1:sum(is.sig.sex))
export(test.sex, "Differential 5-hmC peaks_2-4-10M_Sex_ALLSamples_FDR<0.1_052217.bed")
head(read.table("Differential 5-hmC peaks_2-4-10M_Sex_ALLSamples_FDR<0.1_052217.bed"))

is.sig.age.exp <- tabcom_age_exp$FDR <= 0.10 # No results; doesn't work! (4/10/17 AND 5/22/17)
test.age.exp <- merged$region[is.sig.age.exp]
test.age.exp$score <- -10*log10(tabcom_age_exp$FDR[is.sig.age.exp])
names(test.age.exp) <- paste0("region", 1:sum(is.sig.age.exp))
export(test.age.exp, "Differential 5-hmC peaks_2-4-10M_Age:Exp_ALLSamples_FDR<0.1_052217.bed")
head(read.table("Differential 5-hmC peaks_2-4-10M_Age:Exp_ALLSamples_FDR<0.1_052217.bed"))

#################End of 5/22/17 Overnight Code ################

#######################################################################################
## 5/31/17 - Visualization of genomic coverage facilitated by the extractReads function

#Four regions of exact overlap -- 2-4-10M data(5/30/17):
#chr12 - Bcl11b gene region (intron) overlap between DMR (107965888-107989464) and DHMRs (107972151-107972250;107974201-107974400;107975101-107975300)
#chr13 - Gcnt2 gene region (genes_1to5kb) overlap between DMR (40854768-40858767) and DHMR (40857501-40857750)
#chr17 - Slc22a3 gene region (intron) overlap between DMR (12484898-12489206) and DHMR (12484801-12485100)
#chr19 - Macrod1 gene region (intron) overlap between DMR (7136260-7143741) and DHMRs (7138501-7138700;7141101-7141350)

#Four regions of exact overlap -- ONLY 2M-10M data(4/11/17):
#chr4 - Rgs3 gene region (promoter/exon/intron) overlaps between DMR (62614801, 62635286) and DHMR (62633551, 62633700)
#chr7 - Art1 gene region (intron) overlaps between DMR (102106913, 102114205) and DHMR (102109151, 102109300)
#chr10 - Hmga2 gene region (exon/intron) overlaps between DMR (120373951, 120375570) and DHMR (120375551, 120375700)
#chr10 - Nfic gene region (1to5kb) overlaps between DMR (81432668, 81433772) and DHMR (81433551, 81433800)

#5/31/17 -- Visualize four regions of overlap from both 4/11/17 AND 5/31/17 in GViz: 

#5/30/17:
#chr12: 107972151,107972250; 107974201,107974400; 107975101,107975300
#chr13: 40857501,40857750
#chr17: 12484801,12485100
#chr19: 7138501,7138700; 7141101,7141350

#4/11/17: 
#chr4: 62633551, 62633700
#chr7: 102109151, 102109300
#chr10: 120375551, 120375700
#chr10: 81433551, 81433800

#Note: There was also a single exposure-related region of gene id overlap b/w DMR and DHMR -- 2-4-10M data (5/30/17):
#chr2 - Gnas gene region (intron) DHMR (174315451-174315750) and DMR (174284993-174286537)
#chr2: 174315451,174315750

########################## Esr1 DHM Visualization ########################################
#5/31/17 - Esr1 Region (Matches Mouse Blood Pyro Results)
cur.region <- GRanges("chr10", IRanges(4712147,4712203)) #Esr1 by age
extractReads(bam.files[1], cur.region) 
#Note: No obvious differences in 5-hmC with age. Levels are quite low overall.

#5/31/17 - Esr1 Region 0 (Matches Mouse Blood Pyro Results) #Corresponds to region from tail paper (pyro assay) 
cur.region <- GRanges("chr10", IRanges(4712147,4712203)) #Esr1 by age
extractReads(bam.files[1], cur.region) 
#Note: No obvious differences in 5-hmC with age. Levels are quite low overall.

#8/8/17 - Esr1 Region 1 - Large region that covers a number of sites of interest (See note below)
cur.region <- GRanges("chr10", IRanges(4712147,4712800)) 
extractReads(bam.files[1], cur.region) #Covers age-related promoter and exon DMCs AND three BPA-related DMCs
#5-hmC levels are low in the region covering the DMC at 4712221, but increase in the 
#region covering the DMC at 4712716.
## The full function call to create the GeneRegionTrack from the UCSC data looks like this:
from <- 4712147
to <- 4712800
knownGenes <- UcscTrack(genome = "mm10", chromosome = "chr10",
                        track = "knownGene", from = from, to = to, trackType = "GeneRegionTrack",
                        rstarts = "exonStarts", rends = "exonEnds", gene = "name",
                        symbol = "name", transcript = "name", strand = "strand",
                        fill = "#8282d2", name = "UCSC Genes")
#Furthermore, CpG island annotation can be created using the following:
cpgIslands <- UcscTrack(genome = "mm10", chromosome = "chr10",
                        track = "cpgIslandExt", from = from, to = to,
                        trackType = "AnnotationTrack", start = "chromStart",
                        end = "chromEnd", id = "name", shape = "box",
                        fill = "#006400", name = "CpG Islands")
#Region 1 covers the promoter/exon boundary
pdf("Esr1_DMC_Region1_5hmC_2-4-10M_080817.pdf", height=10, width=10) #8/8/17: Plot created for Esr1 5-hmC
#Need to scale up height and width to get plot to work as pdf
plotTracks(c(itrack, gax, knownGenes, cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() #No obvious differences by exposure


#8/8/17 - Esr1 Region 2
cur.region <- GRanges("chr10", IRanges(4709158,4710157)) #Covers age-related promoter DMC
extractReads(bam.files[1], cur.region) #Very little coverage; not worth investigating the DMC (4709922)

#8/8/17 - Esr1 Region 3 (Figure S4 site from original (Carl Smith Award) draft of paper)
cur.region <- GRanges("chr10", IRanges(4621200,4621600))
extractReads(bam.files[1], cur.region)  #Low-to-moderate coverage; no obvious age or BPA effects

#8/8/17 - Esr1 Region 4 -- Other age:BPA-related DMC (4817892)
cur.region <- GRanges("chr10", IRanges(4817800,4818000))
extractReads(bam.files[1], cur.region) #practically zero coverage. Uninteresting.

#8/8/17 - Esr1 Region 5 - Intronic age- and BPA-related DMCs (4971465, 4972470)
cur.region <- GRanges("chr10", IRanges(4971400,4972500))
extractReads(bam.files[1], cur.region) #Clear 5-hmC peak early on in range (~4971500-4971600)
#Examine region for CpG island annotation and exon/intron boundaries.
from <- 4971400
to <- 4972500
knownGenes <- UcscTrack(genome = "mm10", chromosome = "chr10",
                        track = "knownGene", from = from, to = to, trackType = "GeneRegionTrack",
                        rstarts = "exonStarts", rends = "exonEnds", gene = "name",
                        symbol = "name", transcript = "name", strand = "strand",
                        fill = "#8282d2", name = "UCSC Genes")
#Furthermore, CpG island annotation can be created using the following:
cpgIslands <- UcscTrack(genome = "mm10", chromosome = "chr10",
                        track = "cpgIslandExt", from = from, to = to,
                        trackType = "AnnotationTrack", start = "chromStart",
                        end = "chromEnd", id = "name", shape = "box",
                        fill = "#006400", name = "CpG Islands")
#Region 5 is intronic, but has a clear peak in 5-hmC levels; does not overlap CpG island. Curious!
pdf("Esr1_DMC_Region5_5hmC_2-4-10M_080817.pdf", height=10, width=10) #8/8/17: Plot created for Esr1 5-hmC
#Need to scale up height and width to get plot to work as pdf
plotTracks(c(itrack, gax, knownGenes, cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() #No obvious differences by exposure

#8/8/17 - Esr1 Region 6 - Large region that covers ALL age-, BPA-, and age:BPA-related DMCs
#Also covers the single Exposure-related Esr1 DHMR -- 4694451-4694650
cur.region <- GRanges("chr10", IRanges(4666850,4972500))
extractReads(bam.files[1], cur.region) #Clear 5-hmC peak early on in range (~4971500-4971600)
from <- 4666850
to <- 4972500
knownGenes <- UcscTrack(genome = "mm10", chromosome = "chr10",
                        track = "knownGene", from = from, to = to, trackType = "GeneRegionTrack",
                        rstarts = "exonStarts", rends = "exonEnds", gene = "name",
                        symbol = "name", transcript = "name", strand = "strand",
                        fill = "#8282d2", name = "UCSC Genes")
#Furthermore, CpG island annotation can be created using the following:
cpgIslands <- UcscTrack(genome = "mm10", chromosome = "chr10",
                        track = "cpgIslandExt", from = from, to = to,
                        trackType = "AnnotationTrack", start = "chromStart",
                        end = "chromEnd", id = "name", shape = "box",
                        fill = "#006400", name = "CpG Islands")
#Region 1 covers the promoter/exon boundary
pdf("Esr1_LargeScale_Region6_5hmC_2-4-10M_080817.pdf", height=10, width=10) #8/8/17: Plot created for Esr1 5-hmC
#Need to scale up height and width to get plot to work as pdf
plotTracks(c(itrack, gax, knownGenes, cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() #Although no stunning differences by age or BPA exposure, 5-hmC
#does appear tightly regulated across individuals. VERY INTERESTING!

###################### Gnas DHM Visualization #################################
#5/31/17 - Gnas Region (Exposure-related DHMR; imprinted locus)
cur.region <- GRanges("chr2", IRanges(174315451,174315750)) #Gnas DHMR by exposure
extractReads(bam.files[1], cur.region) 

#6/1/17 - Gnas Region (Exposure-related DMR; imprinted locus)
cur.region <- GRanges("chr2", IRanges(174284993,174286537)) #Gnas by exposure DMR
extractReads(bam.files[1], cur.region) 
#Note: last ~300 bp have more 5-hmC, but it does not seem to vary by exposure; the rest has very little 5-hmC.

#8/8/17 - Zoomed out version of Gnas Region 1 (Exposure-related DHMR; imprinted locus)
cur.region <- GRanges("chr2", IRanges(174315000,174316000)) #Gnas by exposure DHMR
extractReads(bam.files[1], cur.region) 
axistrack <- GenomeAxisTrack(col="black",fontcolor="black",col.border.title="black",range = IRanges(start = 174315451, end = 174315750, names = rep("DHMR")))
plotTracks(axistrack, from = 174315000, to = 174316000, showId = TRUE, col.id="black")
from <- 174315000
to <- 174316000
knownGenes <- UcscTrack(genome = "mm10", chromosome = "chr2",
                        track = "knownGene", from = from, to = to, trackType = "GeneRegionTrack",
                        rstarts = "exonStarts", rends = "exonEnds", gene = "name",
                        symbol = "name", transcript = "name", strand = "strand",
                        fill = "#8282d2", name = "UCSC Genes")
#Furthermore, CpG island annotation can be created using the following:
cpgIslands <- UcscTrack(genome = "mm10", chromosome = "chr2",
                        track = "cpgIslandExt", from = from, to = to,
                        trackType = "AnnotationTrack", start = "chromStart",
                        end = "chromEnd", id = "name", shape = "box",
                        fill = "#006400", name = "CpG Islands")
#1000 bp visualized
pdf("Gnas_DHMR_Region1_5hmC_2-4-10M_Exposure_1000bp_080817.pdf", height=10, width=10) #8/8/17: Plot created for Gnas 5-hmC by age and exposure
plotTracks(c(itrack, axistrack, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() #No obvious differences by exposure outside pre-identified DHMR.


#8/8/17 - Even MORE zoomed out version of Gnas Region 1 (Exposure-related DHMR; imprinted locus)
cur.region <- GRanges("chr2", IRanges(174314000,174317000)) #Gnas by exposure DHMR
extractReads(bam.files[1], cur.region) 
axistrack <- GenomeAxisTrack(col="black",fontcolor="black",col.border.title="black",range = IRanges(start = 174315451, end = 174315750, names = rep("DHMR")))
from <- 174314000
to <- 174317000
knownGenes <- UcscTrack(genome = "mm10", chromosome = "chr2",
                        track = "knownGene", from = from, to = to, trackType = "GeneRegionTrack",
                        rstarts = "exonStarts", rends = "exonEnds", gene = "name",
                        symbol = "name", transcript = "name", strand = "strand",
                        fill = "#8282d2", name = "UCSC Genes")
#Furthermore, CpG island annotation can be created using the following:
cpgIslands <- UcscTrack(genome = "mm10", chromosome = "chr2",
                        track = "cpgIslandExt", from = from, to = to,
                        trackType = "AnnotationTrack", start = "chromStart",
                        end = "chromEnd", id = "name", shape = "box",
                        fill = "#006400", name = "CpG Islands")
#3000 bp visualized
pdf("Gnas_DHMR_Region1_5hmC_2-4-10M_Exposure_3000bp_080817.pdf", height=10, width=10) #8/8/17: Plot created for Gnas 5-hmC by age and exposure
plotTracks(c(itrack, axistrack, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() 

#8/8/17 - The MOST zoomed out version of Gnas -- ENTIRE GENE REGION! 
cur.region <- GRanges("chr2", IRanges(174284306,174346744)) #Gnas gene location from gene card (NCBI) - https://www.ncbi.nlm.nih.gov/gene/14683
extractReads(bam.files[1], cur.region) 
axistrack <- GenomeAxisTrack(col="black",fontcolor="black",col.border.title="black",range = IRanges(start = 174315451, end = 174315750, names = rep("DHMR")))
from <- 174284306
to <- 174346744
knownGenes <- UcscTrack(genome = "mm10", chromosome = "chr2",
                        track = "knownGene", from = from, to = to, trackType = "GeneRegionTrack",
                        rstarts = "exonStarts", rends = "exonEnds", gene = "name",
                        symbol = "name", transcript = "name", strand = "strand",
                        fill = "#8282d2", name = "UCSC Genes")
#Furthermore, CpG island annotation can be created using the following:
cpgIslands <- UcscTrack(genome = "mm10", chromosome = "chr2",
                        track = "cpgIslandExt", from = from, to = to,
                        trackType = "AnnotationTrack", start = "chromStart",
                        end = "chromEnd", id = "name", shape = "box",
                        fill = "#006400", name = "CpG Islands")
#Entire gene (~162,000 bp) visualized
pdf("Gnas_Entire_Gene_5hmC_2-4-10M_080817.pdf", height=10, width=10) #8/8/17: Plot created for Gnas 5-hmC by age and exposure
plotTracks(c(itrack, axistrack, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() 

#20 bp window version of GNAS data (1/8/18):
pdf("Gnas_Entire_Gene_5hmC_20bpWindow_2-4-10M_010818.pdf", height=10, width=10) #8/8/17: Plot created for Gnas 5-hmC by age and exposure
plotTracks(c(itrack, axistrack, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() 
#Note: Looked very similar to 100 bp window size data; peaks were better defined, but in the same places.

################### Igf2/H19 gene region visualized #####################
#For curiosity's sake, plot the Igf2 region:

##Plots are made using the Gviz package (Note: Need "data" defined in step 1 to proceed)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)


#11/1/17 - The MOST zoomed out version of Igf2 -- ENTIRE GENE REGION! 
cur.region <- GRanges("chr7", IRanges(142650768,142658804)) #Igf2 gene location from UCSC genome browser (mm10); size = 8307 bp (including UTRs)
extractReads(bam.files[1], cur.region) 
from <- 142650768
to <- 142658804

#11/1/17 - Zoom in on potential region of differential 5-hmC by BPA (chr7: 142652500,142654200)
cur.region <- GRanges("chr7", IRanges(142652500,142654200)) #Igf2 gene location from UCSC genome browser (mm10); size = 8307 bp (including UTRs)
extractReads(bam.files[1], cur.region) 
from <- 142652500
to <- 142654200

#11/1/17 - Zoom in even further on potential region of differential 5-hmC by BPA (chr7: 142653200,142654200)
cur.region <- GRanges("chr7", IRanges(142653450,142654200)) #Igf2 gene location from UCSC genome browser (mm10); size = 8307 bp (including UTRs)
extractReads(bam.files[1], cur.region) 
from <- 142653450
to <- 142654200

#11/1/17 - The MOST zoomed out version of H19 -- ENTIRE GENE REGION! 
cur.region <- GRanges("chr7", IRanges(142575530,142578146)) #H19 gene location from UCSC genome browser (mm10); size = 2617 bp
extractReads(bam.files[1], cur.region) 
from <- 142575530
to <- 142578146

#11/1/17 - Zoom in on potential H19 region of differential 5-hmC by BPA (chr7: 142577250,142577750)
cur.region <- GRanges("chr7", IRanges(142577250,142577750)) #H19 gene location from UCSC genome browser (mm10); size = 2617 bp
extractReads(bam.files[1], cur.region) 
from <- 142577250
to <- 142577750

#11/1/17 - Plagl1 imprinted region -- ENTIRE GENE (The DHMR seems cut-off on both ends)
#Exposure-related DHMR (exonic)
#chr10: 13090788-13131695
cur.region <- GRanges("chr10", IRanges(13090788,13131695)) #Plagl1 gene location from UCSC genome browser (mm10); size = 40908 bp 
extractReads(bam.files[1], cur.region) 
from <- 13090788
to <- 13131695


knownGenes <- UcscTrack(genome = "mm10", chromosome = "chr10",
                        track = "knownGene", from = from, to = to, trackType = "GeneRegionTrack",
                        rstarts = "exonStarts", rends = "exonEnds", gene = "name",
                        symbol = "name", transcript = "name", strand = "strand",
                        fill = "#8282d2", name = "UCSC Genes")
#Furthermore, CpG island annotation can be created using the following:
cpgIslands <- UcscTrack(genome = "mm10", chromosome = "chr10",
                        track = "cpgIslandExt", from = from, to = to,
                        trackType = "AnnotationTrack", start = "chromStart",
                        end = "chromEnd", id = "name", shape = "box",
                        fill = "#006400", name = "CpG Islands")

#Define plotting range and establish reference to HMeDIP-seq data:
chr <- as.character(unique(seqnames(cur.region))) #define chr location
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene #establish mm10 as genome of interest
seqlevels(txdb) #view all levels of mm10 included from package
tx <- transcripts(txdb) #restrict to list
seqlevels(tx, force=TRUE) <- c("chr1") ## Drop all items in list except 'chr1'.
seqlevels(tx) #view levels
gen <- genome(tx) #generate genome from restricted list

itrack <- IdeogramTrack(genome=gen, chr=chr) #Creates ideogram to place at top of the graph to indicate location on chromosome

collected <- list() #Note: generates collection of regions for ploting
for (i in 1:length(bam.files)) {
  reads <- extractReads(bam.files[i], cur.region)
  pcov <- as(coverage(reads[strand(reads)=="+"])/data$totals[i]*1e6, "GRanges")
  ncov <- as(coverage(reads[strand(reads)=="-"])/data$totals[i]*1e6, "GRanges")
  ptrack <- DataTrack(pcov, type="histogram", lwd=0, fill=rgb(0,0,1,.4), ylim=c(0,0.4),
                      name="5-HmC", col.axis="black", col.title="black")
  ntrack <- DataTrack(ncov, type="histogram", lwd=0, fill=rgb(1,0,0,.4), ylim=c(0,0.4))
  collected[[i]] <- OverlayTrack(trackList=list(ptrack,ntrack))
} 

#Entire Igf2 gene + UTRs (~8,000 bp) visualized
pdf("Igf2_Entire_Gene_5hmC_2-4-10M_110117.pdf", height=10, width=10) #11/1/17: Plot created for Igf2 5-hmC by age and exposure
plotTracks(c(itrack, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() 

#Zoom in on potential Igf2 region of differential 5-hmC (chr7: 142652500,142654200)
pdf("Igf2_1700bp_DHMR_5hmC_2-4-10M_110117.pdf", height=10, width=10) #11/1/17: Plot created for Igf2 5-hmC by age and exposure
plotTracks(c(itrack, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() 

#Zoom in even further on potential Igf2 region of differential 5-hmC (chr7: 142653450,142654200)
pdf("Igf2_750bp_DHMR_5hmC_2-4-10M_110117.pdf", height=10, width=10) #11/1/17: Plot created for Igf2 5-hmC by age and exposure
plotTracks(c(itrack, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() 

#Entire H19 gene (~2600 bp) visualized
pdf("H19_Entire_Gene_5hmC_2-4-10M_110117.pdf", height=10, width=10) #11/1/17: Plot created for H19 5-hmC by age and exposure
plotTracks(c(itrack, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() 

#Zoom in on potential H19 gene of differential 5-hmC by BPA (~500 bp) (chr7: 142577250,142577750)
pdf("H19_500bp_5hmC_2-4-10M_110117.pdf", height=10, width=10) #11/1/17: Plot created for H19 5-hmC by age and exposure
plotTracks(c(itrack, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() 

#11/1/17 - Plagl1 imprinted region -- ENTIRE GENE (The DHMR seems cut-off on both ends)
#Exposure-related DHMR (exonic)
#chr10: 13090788-13131695
pdf("Plagl1_Entire_Gene_5hmC_2-4-10M_110117.pdf", height=10, width=10) #11/1/17: Plot created for Plagl1 5-hmC by age and exposure
plotTracks(c(itrack, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))
dev.off() 

#12/4/17 - Phactr2	(phosphatase and actin regulator 2) imprinted region (INCREASED Y LIMIT)
#Exposure-related DHMR (intronic)
#chr10: 13402051,13402250
cur.region <- GRanges("chr10", IRanges(13402051,13402250)) #Phactr2 gene location from UCSC genome browser (mm10); size = 40908 bp 
extractReads(bam.files[1], cur.region) 
gax <- GenomeAxisTrack(col="black")
from <- 13402051
to <- 13402250
#Create pdf of visualization:
pdf("Phactr2_5hmC_2-4-10M_Exposure_ALLSamples_120417.pdf", height=10, width=10) #12/4/17: Plot created for Phactr2 5-hmC by age and exposure
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off() 

#Note: adjusted ylim up to 120% (1.2) for this region's visualization - 04/03/17
#Note: adjusted ylim down to 40% (0.4) for region visualization - 04/04/17
#Note: adjusted ylim down to 30% (0.3) for region visualization - 04/12/17
gax <- GenomeAxisTrack(col="black")
#Note: chr and gen defined from cur.region above!
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
plotTracks(c(itrack, gax, knownGenes, collected), from=start(cur.region), to=end(cur.region))
plotTracks(c(itrack, gax, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))

#################### Imprinted gene regions with DHMRs #############################
#6/15/17 - Igf1 Gene Exposure-related DHMR (intronic):
#chr10: 87929451,87930200

#6/15/17 - Ppp1r9a (protein phosphatase 1, regulatory (inhibitor) subunit 9A) imprinted region
#Exposure-related DHMR (intronic)
#chr6:5079501,5079850

#6/15/17 - Plagl1 (pleiomorphic adenoma gene-like 1) imprinted region
#Exposure-related DHMR (exonic)
#chr10: 13130651, 13131050

#11/1/17 - Plagl1 imprinted region -- ENTIRE GENE (The DHMR seems cut-off on both ends)
#Exposure-related DHMR (exonic)
#chr10: 13090788, 13131695


#6/15/17 - Phactr2	(phosphatase and actin regulator 2) imprinted region
#Exposure-related DHMR (intronic)
#chr10: 13402051,13402250

#12/4/17 - Phactr2	(phosphatase and actin regulator 2) imprinted region (INCREASED Y LIMIT)
#Exposure-related DHMR (intronic)
#chr10: 13402051,13402250

#6/15/17 - Pde4d (phosphodiesterase 4D, cAMP specific)
#Exposure-related DHMR (intronic)
#chr13: 108963951,108964200

#6/15/17 - Pde10a (phosphodiesterase 10A)
#Exposure-related DHMR (intronic)
#chr17: 8746901,8747050

#6/15/17 - Klf14(Kruppel-like factor 14)
#Exposure-related DHMR (exonic)
#chr6: 30957751,30958050

#6/15/17 - Kcnq1 (potassium voltage-gated channel, subfamily Q, member 1)
#Exposure-related DHMR (intronic)
#chr6: 30957751,30958050

#6/15/17 - 	Grb10 (growth factor receptor bound protein 10)
#Exposure-related DHMR (intronic)
#chr11:12004801,12005100

#6/15/17 - 	Cmah (cytidine monophospho-N-acetylneuraminic acid hydroxylase)
#Exposure-related DHMR (intronic)
#chr13:24335001,24335150

#6/15/17 - 	Airn (antisense Igf2r RNA)
#Exposure-related DHMR (intronic)
#chr17:12828001,12828550

#6/15/17 -- Snrpn (Imprinted locus) Exposure-related DHMR (intron):
#chr7: 60207401,60207600

#5/30/17 -- four regions of DMR/DHMR overlap in GViz: 
#Bcl11b -- #chr12: 107972151,107972250
cur.region <- GRanges("chr12", IRanges(107972151,107972250)) #Bcl11b by age
extractReads(bam.files[1], cur.region)
#Gcnt2 -- chr13: 40857501,40857750
cur.region <- GRanges("chr13", IRanges(40857501,40857750)) #Gcnt2 by age
extractReads(bam.files[1], cur.region)
#Slc22a3 -- chr17: 12484801,12485100
cur.region <- GRanges("chr17", IRanges(12484801,12485100)) #Slc22a3 by age
extractReads(bam.files[1], cur.region)
#Macrod1 -- chr19: 7138501,7138700
cur.region <- GRanges("chr19", IRanges(7138501,7138700)) #Macrod1 by age
extractReads(bam.files[1], cur.region)

#4/11/17 - four regions of DHMR and DMR overlap in GViz:
#Rgs3 -- chr4: 62633551, 62633700
cur.region <- GRanges("chr4", IRanges(62633551,62633700)) #Rgs3 by age
extractReads(bam.files[1], cur.region)
#Art1 -- chr7: 102109151, 102109300
cur.region <- GRanges("chr7", IRanges(102109151,102109300)) #Art1 by age
extractReads(bam.files[1], cur.region)
#Hmga2 -- chr10: 120375551, 120375700
cur.region <- GRanges("chr10", IRanges(120375551,120375700)) #Hmga2 by age
extractReads(bam.files[1], cur.region)
#Nfic -- chr10: 81433551, 81433800
cur.region <- GRanges("chr10", IRanges(81433551,81433800)) #Nfic by age
extractReads(bam.files[1], cur.region)

##Plots are made using the Gviz package (Note: Need "data" defined in step 1 to proceed)
library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

chr <- as.character(unique(seqnames(cur.region))) #define chr location
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene #establish mm10 as genome of interest
seqlevels(txdb) #view all levels of mm10 included from package
tx <- transcripts(txdb) #restrict to list
seqlevels(tx, force=TRUE) <- c("chr1") ## Drop all items in list except 'chr1'.
seqlevels(tx) #view levels
gen <- genome(tx) #generate genome from restricted list

itrack <- IdeogramTrack(genome=gen, chr=chr) #Creates ideogram to place at top of the graph to indicate location on chromosome

collected <- list() #Note: generates collection of regions for ploting
for (i in 1:length(bam.files)) {
  reads <- extractReads(bam.files[i], cur.region)
  pcov <- as(coverage(reads[strand(reads)=="+"])/data$totals[i]*1e6, "GRanges")
  ncov <- as(coverage(reads[strand(reads)=="-"])/data$totals[i]*1e6, "GRanges")
  ptrack <- DataTrack(pcov, type="histogram", lwd=0, fill=rgb(0,0,1,.4), ylim=c(0,0.3),
                      name="5-HmC", col.axis="black", col.title="black")
  ntrack <- DataTrack(ncov, type="histogram", lwd=0, fill=rgb(1,0,0,.4), ylim=c(0,0.3))
  collected[[i]] <- OverlayTrack(trackList=list(ptrack,ntrack))
} 
#Note: adjusted ylim up to 120% (1.2) for this region's visualization - 04/03/17
#Note: adjusted ylim down to 40% (0.4) for region visualization - 04/04/17
#Note: adjusted ylim down to 30% (0.3) for region visualization - 04/12/17
gax <- GenomeAxisTrack(col="black")
#Note: chr and gen defined from cur.region above!
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
plotTracks(c(itrack, gax, knownGenes, collected), from=start(cur.region), to=end(cur.region))
plotTracks(c(itrack, gax, knownGenes,cpgIslands, collected), from=start(cur.region), to=end(cur.region))


#PDF version of Esr1 Region DHMR R Plot to ONES Server (for manipulation in Illustrator):
pdf("Esr1_5hmC_2-4-10M_Age_ALL Samples_053117.pdf") #5/31/17: Plot created for Esr1 5-hmC by age
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off() #No apparent shifts in Esr1 5-hmC with age or exposure.

#PDF version of Gnas Region DHMR R Plot to ONES Server (for manipulation in Illustrator):
pdf("Gnas_5hmC_2-4-10M_Age_ALL Samples_053117.pdf") #5/31/17: Plot created for Gnas 5-hmC by age
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off() #WOW! HUGE change in 5-hmC coverage based on BPA exposure!!!!

pdf("Gnas_DMR Region_5hmC_2-4-10M_Age_ALL Samples_060117.pdf") #5/31/17: Plot created for Gnas 5-hmC by age
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off() #No obvious differences by exposure;however, the last few 100 bps look interesting (zoom in further!)



###############2-4-10M Regions of DHMR/DMR overlap (5/30/17):
#Save Age-based Bcl11b Region DHMR R plot to ONES server:
pdf("Bcl11b_5hmC_2-4-10M_Age_ALL Samples_053117.pdf") #5/31/17: Plot created for Bcl11b 5-hmC by age
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off() #Age-related results are subtle; not visually identifiable

#Save Age-based Gcnt2 Region DHMR R plot to ONES server:
pdf("Gcnt2_5hmC_2-4-10M_Age_ALL Samples_053117.pdf") #5/31/17: Plot created for Gcnt2 5-hmC by age
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off() #Age-related results are subtle; not visually identifiable

#Save Age-based Slc22a3 Region DHMR R plot to ONES server:
pdf("Slc22a3_5hmC_2-4-10M_Age_ALL Samples_053117.pdf") #5/31/17: Plot created for Slc22a3 5-hmC by age
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off() #Age-related results are subtle; not visually identifiable

#Save Age-based Macrod1 Region DHMR R plot to ONES server:
pdf("Macrod1_5hmC_2-4-10M_Age_ALL Samples_053117.pdf") #5/31/17: Plot created for Macrod1 5-hmC by age
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off() #Age-related results are subtle; not visually identifiable

###############2M-10M Regions of DHMR/DMR overlap (4/11/17):
#Save Age-based Rgs3 Region DHMR R plot to ONES server:
pdf("Rgs3_5hmC_2-4-10M_Age_ALL Samples_053117.pdf") #5/31/17: Plot created for Rgs3 5-hmC by age
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off()

#Save Age-based Art1 Region DHMR R plot to ONES server:
pdf("Art1_5hmC_2-4-10M_Age_ALL Samples_053117.pdf") #4/18/17: Plot created for Art1 5-hmC by age
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off()

#Save Age-based Hmga2 Region DHMR R plot to ONES server:
pdf("Hmga2_5hmC_2-4-10M_Age_ALL Samples_053117.pdf") #4/18/17: Plot created for Hmga2 5-hmC by age
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off()

#Save Age-based Nfic Region DHMR R plot to ONES server:
pdf("Nfic_5hmC_2-4-10M_Age_ALL Samples_053117.pdf") #4/18/17: Plot created for Nfic 5-hmC by age
plotTracks(c(itrack, gax, collected), from=start(cur.region), to=end(cur.region))
dev.off()

################## End of 5/31/17 5-hmC Visualization ####################
