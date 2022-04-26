################################################
##   Bioinf545/ Biostats646/ Stats545
##   Feb. 26, 2016
##   Lab 4: ChIP-seq analysis & clustering
################################################


head /class/data/bio545w16/lab4/chr19_Lane1.bed
macs2 callpeak
cd ~ 

### Run MACS peak-finder on ATF4 data from chromosome 19
##  More commonly for a chip-seq dataset, you'll be starting with BAM input rather than BED.
##  -t specifies the ChIP sample;  -c specifies the control sample
###  P-value cutoff of 10^(-3)
macs2 callpeak -t /class/data/bio545w16/lab4/chr19_Lane1.bed -c /class/data/bio545w16/lab4/chr19_Lane2.bed --format=BED  --gsize=61342430  --pvalue=1e-3  --name=ATF4e-3

### q-value cutoff of 0.05
macs2 callpeak -t /class/data/bio545w16/lab4/chr19_Lane1.bed -c /class/data/bio545w16/lab4/chr19_Lane2.bed -f=BED --gsize=61342430 --qvalue=.05 --name=ATF4q05

###  q-value cutoff of 0.01
macs2 callpeak -t /class/data/bio545w16/lab4/chr19_Lane1.bed -c /class/data/bio545w16/lab4/chr19_Lane2.bed -f=BED --gsize=61342430 --qvalue=.01 --name=ATF4q01

head /class/data/bio545w16/lab4/ATF4q01_peaks.narrowPeak

### View top rows of the MACS output bed file for all of the ATF4 data
head /class/data/bio545w16/lab4/region_ATF4.bed


####################################################################################
### Start RStudio or R on your computer to read in the ChIP-seq results and do 
###  some post-processing. 
####################################################################################

library(ChIPpeakAnno)   # Load the ChIP peak annotation package
chopbed<-read.delim("region_CHOP.bed", header=F, stringsAsFactors=F) # Read in peaks for whole genome
atf4bed<-read.delim("region_ATF4.bed", header=F, stringsAsFactors=F)
#chopbed<-read.delim("/class/data/bio545w16/lab4/region_CHOP.bed", header=F, stringsAsFactors=F)
#atf4bed<-read.delim("/class/data/bio545w16/lab4/region_ATF4.bed", header=F, stringsAsFactors=F)
atf4bed[1:4,]
dim(chopbed)  # Show number of rows = how many peaks
dim(atf4bed)

###  Read in the results from MACS for chr19 (using q-value < 0.01 cutoff)
atf4.macs.c19<-read.delim("ATF4q01_peaks.narrowPeak", header=F, stringsAsFactors=F)
#atf4.macs.c19<-read.delim("/class/data/bio545w16/lab4/ATF4q01_peaks.narrowPeak", header=F, stringsAsFactors=F)
# Note, for some files you would need to use skip=1 (or some other number) if there's 1 or more extra lines of information 
# at the top of the file.
dim(atf4.macs.c19)  # should be 401 peaks

colnames(atf4.macs.c19)<-c('chr','start','stop','peak','score','na','fc','nlog10p','nlog10q','summit')
atf4.macs.c19[1:4,]

###  MACS negative log10 q-values
summary(atf4.macs.c19$nlog10q)  # negative log base 10 q-values
10^(-3.202)  # Transform to get back to actual q-values
10^(-254.2)

##  Compare number of peaks on chr19 (which peak caller called more peaks?)
atf4.c19<-atf4bed[atf4bed[,1]=="chr19",]  # Extract out just the peaks on chromosome 19
dim(atf4.c19)

##  Compare peak lengths between MACS and Erange
macs.lengths<-atf4.macs.c19[,3]-atf4.macs.c19[,2]  # Calculate width of peaks called by MACS
E.lengths<-atf4.c19[,3]-atf4.c19[,2]  # Calculate width of peaks called by Erange
summary(macs.lengths)
summary(E.lengths)
#png("peak_length_histograms.png",width=500,height=300)
par(mfrow=c(1,2))
hist(macs.lengths)
hist(E.lengths)
#dev.off()

###  Examine the fold changes
summary(atf4.macs.c19$fc)
#png("atf4.macs.c19 fc histogram.png",width=300,height=300)
hist(atf4.macs.c19$fc,breaks=32)
#dev.off()

#png("atf4-macs-FC_vs_nlog10p.png",width=300,height=300)
plot(atf4.macs.c19$fc,atf4.macs.c19$nlog10p,xlab="fold change",ylab="-log10(p-value)")
#dev.off()

# Using IRanges and Ranged data to add annotation and detect overlapping peaks
ir <- IRanges(c(103,200,218,311), width=c(10,21,13,9),names=c('a','b','c','d'))
ir
start(ir)  # the start sites
end(ir)	# the end sites
rev(ir) # reverse order
ir[[1]]  # get base locations for first site

chopr<-BED2RangedData(chopbed)
atf4r<-BED2RangedData(atf4bed)
chopr[1:4,]
chop.spaces<-IRanges::space(chopr)  # Force R to use the space() function from the IRanges package. Gives the chromosome of each peak

table(chop.spaces) # How many peaks on each chromosome?
length(chopr)  # the number of chromosomes with at least one peak
dim(chopr)   # total number of peaks across genome

myexon=c(0,1,0,1) # in exon?  
myspace=c('chr1','chr2', 'chr2', 'chr1')
myGC=runif(4)  # proportion of GC’s in each range
myexon
myspace
myGC
rd=RangedData(ranges=ir, space=myspace, GC=myGC, exon=myexon)
?RangedData
  
rd  # Notice R sorted them automatically
start(rd) 	# Get start positions
names(rd) 	# Get chromosome names
ranges(rd) 	# Get the ranges part of object
values(rd)	# Get GC and exon info (value columns)
rdd<-as.data.frame(rd) # convert to a data frame
rdd
reduce(rd) # Merges any overlapping ranges; now there should only be 3 ranges.


# Getting back to our ChIP-seq data... findOverlaps between peak callers for chr19
c19.erange<-BED2RangedData(atf4.c19)
c19.macs<-BED2RangedData(atf4.macs.c19[,1:5])
e.macs19<-findOverlaps(c19.erange,c19.macs)
e.macs19<-as.matrix(e.macs19)
dim(e.macs19)  # 104
e.macs19[1:5,]
atf4.c19[e.macs19[1:5,1],]
atf4.macs.c19[e.macs19[1:5,2],]


# Annotate peaks to transcription start sites (TSSs)
data(package='ChIPpeakAnno')
data(TSS.mouse.NCBIM37)  # This line may take a few minutes to run...
annotPeak<-annotatePeakInBatch(chopr[1:6,], AnnotationData=TSS.mouse.NCBIM37)
peak6<-as.data.frame(annotPeak)
peak6

# Annotate peaks to microRNAs
mart = useMart(biomart="ensembl", dataset= "mmusculus_gene_ensembl")
miRAnnot = getAnnotation(mart, featureType="miRNA")
annotmiR.peak = annotatePeakInBatch(chopr, AnnotationData=miRAnnot)
miR.peak<-as.data.frame(annotmiR.peak)
miR.peak.sorted<-miR.peak[order(miR.peak$shortestDistance),]
miR.peak.sorted[1:5,]  # What are the closest peaks to a miRNA?

# Find overlapping peaks using the ChIPpeakAnno package
overlaps<-findOverlapsOfPeaks(chopr,atf4r, maxgap=100)
names(overlaps)  # Show what is in this object
attributes(overlaps$overlappingPeaks)  # Show the contents of this object
chopr.atf4r<-overlaps$overlappingPeaks$"chopr///atf4r"
chopr.atf4r[1:3,]  #  Top of table of annotated results
dim(chopr.atf4r)

?findOverlapsOfPeaks  # Check out other options

makeVennDiagram(RangedDataList(chopr,atf4r), NameOfPeaks=c("CHOP","ATF4"), maxgap=0, totalTest=300000, cex=2)
fisher.test(matrix(c(1425,1173,1598,295804), 2,2))
fisher.test(matrix(c(1425,1173,1598,50000), 2,2))  # Is it still significant if we assume a much smaller total number of windows tested?

