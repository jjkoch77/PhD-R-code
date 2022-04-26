##############################################################################
##  Bioinf545 / Biostat646 / Stats545
##  Lab 4 ;  2-26-16
##  This code file shows how to use the R package 'annotatr'
##  to annotate the peaks in a bed file to different types of genomic regions
##
##  The vignette for annotatr is at: https://github.com/rcavalcante/annotatr (scroll down)
##  Working through this file is optional; annotatr is still in "beta version".
##  However, you might find it more useful than ChIPpeakAnno.
##  Please report any bugs to Raymond at rcavalca@umich.edu
##  Send questions to Maureen- sartorma@umich.edu
###############################################################################

# Working directory
getwd()
setwd('~/Desktop') # Change to wherever you saved the file ATF4q01_peaks.narrowPeak.

##  If you don't have the devtools package installed, install it now.
# Load devtools, which allows easy installation of an R package from github
library(devtools)

# Install dependencies from Bioconductor (devtools doesn't support auto-install of these)
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocStyle","GenomeInfoDb","IRanges","GenomicRanges"))

# Install annotatr and load it
devtools::install_github('rcavalcante/annotatr')
library(annotatr)

# List the supported annotations
supported_annotations()
## The following code uses the detailed genes annotation for mouse, and shows you how to use a custom annotation.
## Go to https://github.com/rcavalcante/annotatr  and scroll down to see an illustration of 
##    the "UCSC knownGenes" detailedgenes annotation.  
##    "UCSC knowngenes" is a commonly used gene definition file.


##  A feature of annotatr is that it allows annotating regions to any custom annotation dataset.
##  We'll show how we can use this feature to annotate to another ChIP-seq peak set 
##  to find overlaps.

##  Go to  http://mouseencode.org > Data Search > Select ChIP-seq + mm9 + bed broadPeak (along the left)
##  This gets you to a list of 205 experiments to choose from.
##  Choose a broadpeak file (e.g. H3K4me3 in limb), and under the "Processed data" section, download the bed broadPeak.
##  Save to your home directory on bcs2

# Use the following bash code on command line on bcs2 to reformat the broadPeak file for read_annotations()
# Replace the file names below with the one that you chose.
gunzip -c ENCFF438EEG.bed.gz | head
# To see what it looks like

cut -f 1-3 <(gunzip -c ENCFF438EEG.bed.gz) > NANOG_mm9.bed
## cut is a UNIX utility that lets you grab columns from tab-delimited file (tab by default, but this can be changed). 
##  The <(gunzip...) part is a subshell that unzips the file to standard out to pipe into cut. 
## Now transfer the NANOG_mm9.bed (or whatever name you called it) back to your computer and save it in 
##   the same directory as ATF4q01_peaks.narrowPeak

# Back in R...  Read in the custom annotation from mouse ENCODE (Change the file name if necessary) 
mm9_custom_ENCODE = read_annotations(
	file = 'NANOG_mm9.bed',
	genome = 'mm9',
	annotation_name = 'ENCODE')

################################
# Analysis for ATF4 q01 (chr19) in mouse

# Read file. It's a good idea to name the columns
atf4gr <- read_bed(file = "ATF4q01_peaks.narrowPeak", genome = "mm9",
	col.names = c('chrom','start','end','name','score','strand','signalValue','pvalue','qvalue','peak'))

# Annotate the atf4 regions to the detailed genes annotation of mm9 and the custom dataset regions
annots <- annotate_regions(
	regions = atf4gr,
	annotations=c('mm9_detailedgenes','mm9_custom_ENCODE'),
	use.score=T)
annots[1:4,]
#  annots is a special type of dataframe using the dplyr package.
# By default, dplyr::tbl_df objects have nice printing properties, but it
# hides extra columns that would ordinarily wrap. You can see them all with:
head(as.data.frame(annots))

## You could write the annotations to a tab-delimited text file with:
write.table(as.data.frame(annots),file="annotatr_example_annotations.txt",sep="\t",quote=F)
## You could then open this file with Excel for further processing.

# Summarize the number of regions per annotation type
counts <- summarize_annotations(annots)
print(counts)

# Set the order of the annotations
annot_order = c('mm9_custom_ENCODE',
				'mm9_knownGenes_1to5kb',
				'mm9_knownGenes_promoters',
				'mm9_knownGenes_exons5UTRs',
				'mm9_knownGenes_introns5UTRs',
				'mm9_knownGenes_exonsCDSs',
				'mm9_knownGenes_intronsCDSs',
				'mm9_knownGenes_exons3UTRs',
				'mm9_knownGenes_introns3UTRs')

# Plot the number of regions per annotation type
plot_counts <- visualize_annotation(
	annotated_regions = annots, annotation_order = annot_order,
	plot_title = 'ATF4 Binding in mm9', x_label = 'Annotations', y_label = '# Regions')
plot_counts

# Plot the number of regions per pair of annotation type
# Note where the ENCODE overlaps are: promoters, exons in 5'UTRs, CDS introns...?
plot_coannotations <- visualize_coannotations(
	annotated_regions = annots, annotation_order = annot_order,
  	plot_title = 'ATF4 Binding in Pairs of Annotations', axes_label = 'Annotations')
plot_coannotations

# Additional plots are described in the vignette.


##############################
# Analysis for ATF4 sites across the genome (from Erange peak finder) in mouse

# Read file. It's a good idea to name the columns
atf4gr <- read_bed(file = "region_ATF4.bed", genome = "mm9",
	col.names = c('chrom','start','end','name'))

# Annotate the atf4 regions to the detailed genes annotation of mm9
annots <- annotate_regions(
	regions = atf4gr,
	annotations=c('mm9_detailedgenes','mm9_custom_ENCODE'),
	use.score=F)

# By default, dplyr::tbl_df objects have nice printing properties, but it
# hides extra columns that would ordinarily wrap. You can see them all with:
head(as.data.frame(annots))

# Summarize the number of regions per annotation type
counts <- summarize_annotations(annots)
print(counts)

# Set the order of the annotations
annot_order = c('mm9_custom_ENCODE',
				'mm9_knownGenes_1to5kb',
				'mm9_knownGenes_promoters',
				'mm9_knownGenes_exons5UTRs',
				'mm9_knownGenes_introns5UTRs',
				'mm9_knownGenes_exonsCDSs',
				'mm9_knownGenes_intronsCDSs',
				'mm9_knownGenes_exons3UTRs',
				'mm9_knownGenes_introns3UTRs')

# Plot the number of regions per annotation type
plot_counts <- visualize_annotation(
	annotated_regions = annots, annotation_order = annot_order,
	plot_title = 'ATF4 Binding in mm9', x_label = 'Annotations', y_label = '# Regions')
plot_counts

# Plot the number of regions per pair of annotation type
# Note which ENCODE binding sites are also in which annotations...
plot_coannotations <- visualize_coannotations(
	annotated_regions = annots, annotation_order = annot_order,
  	plot_title = 'ATF4 Binding in Pairs of Annotations', axes_label = 'Annotations')
plot_coannotations

