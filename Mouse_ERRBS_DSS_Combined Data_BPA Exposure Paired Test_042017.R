##Running DSS to examine paired differential methylation in mouse blood eRRBS data

##Install dependency (BSseq package):
source("https://bioconductor.org/biocLite.R")
biocLite("bsseq")
library(bsseq)
biocLite("DSS")
library(DSS)
biocLite("annotatr")
library(annotatr)
biocLite("csaw")
library(csaw)

##Extract methylation counts using bismark_methylation_extractor function:
#bismark_methylation_extractor -s -bedGraph reads.fastq_bismark.sam. 
#This will create multiple txt
#files to summarize methylation call and cytosine context, a bedGraph file to 
#display methylation percentage, and a coverage file containing counts information. 
#The count file contain following columns:chr, start, end,
#methylation%, count methylated, count unmethylated. 
#This file can be modified to make the input file for DSS.

##Note: coverage files already exist on ONES server. Copy them over, then reformat to
##fit the necessary DSS format:

#DSS requires data from each BS-seq experiment to be summarized into following information for each CG position:
#chromosome number, genomic coordinate, total number of reads, and number of reads showing methylation. For a
#sample, this information are saved in a simple text file, with each row representing a CpG site. Below shows an example
#of a small part of such a file:
#  chr pos N X
#chr18 3014904 26 2
#chr18 3031032 33 12
#chr18 3031044 33 13

##Read in .cov bismark files for each sample
##Example:
#read.bismark(files, sampleNames, rmZeroCov = FALSE, strandCollapse = TRUE, fileType = c("cov", "oldBedGraph", "cytosineReport"), mc.cores = 1, verbose = TRUE)

#files - Input files. Each sample is in a different file. Input files are created by running Bismark's methylation extractor; see Note for details.
#sampleNames - sample names, based on the order of files.
#rmZeroCov - Should methylation loci that have zero coverage in all samples be removed. This will result in a much smaller object if the data originates from (targeted) capture bisulfite sequencing.
#strandCollapse - Should strand-symmetric methylation loci, e.g., CpGs, be collapsed across strands. This option is only available if fileType = "cytosineReport" since the other file types do not contain the necessary strand information.
#fileType - The format of the input file; see Note for details.
#mc.cores - The number of cores used. Note that setting mc.cores to a value greater than 1 is not supported on MS Windows, see the help page for mclapply.
#verbose - Make the function verbose.

#Set working directory
getwd()
setwd("/jjkoch/03-methylation_call")


##################### Covariate experimental design #################

#In DSS, BS-seq data from a general experimental design (such as crossed experiment, or experiment with covariates) is
#modeled through a generalized linear model framework. 
#We use "arcsine"link function instead of the typical logit link for it better deals with data at boundaries (methylation levels close to 0 or 1).
#Linear model fitting is done through ordinary least square on transformed methylation levels. 
#Standard errors for the estimates are derived with consideration of count data distribution and transformation. 
#A Wald test is applied to perform hypothesis testing.

################Step 1##################

##Establish combined 2M, 4M, and 10M dataset for BPA comparison with max sample size - 4/20/17 
bsseq_combined <- read.bismark(files = c("171c_2M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "241c_2M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "411b_2M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "441d_2M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "981c_2M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "991d_2M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "161a_2M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "101a_2M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "501a_2M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "431a_2M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "941a_2M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "951a_2M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "171c_4M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "241c_4M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "411b_4M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "441d_4M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "981c_4M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "991d_4M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "161a_4M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "101a_4M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "501a_4M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "431a_4M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "941a_4M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "951a_4M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "171c_10M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "241c_10M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "411b_10M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "441d_10M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "981c_10M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "991d_10M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "161a_10M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "101a_10M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "501a_10M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "431a_10M_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "941a_10M_R1_trimmed.fq.gz_bismark_bt2.bismark.cov",
                                         "951a_10M_trimmed.fq.gz_bismark_bt2.bismark.cov" ),
                               sampleNames = c("2M_1", "2M_2","2M_3","2M_4","2M_5","2M_6",
                                               "2M_7", "2M_8","2M_9","2M_10","2M_11","2M_12",
                                               "4M_1", "4M_2","4M_3","4M_4","4M_5","4M_6",
                                               "4M_7", "4M_8","4M_9","4M_10","4M_11","4M_12",
                                               "10M_1","10M_2","10M_3","10M_4","10M_5", "10M_6",
                                               "10M_7","10M_8","10M_9","10M_10","10M_11", "10M_12"),
                               rmZeroCov = TRUE,
                               strandCollapse = FALSE,
                               fileType = "cov",
                               verbose = TRUE)

############ Step 2 ###################

##Must establish design prior to running multifactor function.
##all 2M,4,and 10M samples combined in single matrix 
n <- c(1:36)
age <- c("2M", "2M", "2M", "2M", "2M", "2M","2M", "2M", "2M", "2M", "2M", "2M",
         "4M", "4M", "4M", "4M", "4M", "4M","4M", "4M", "4M", "4M", "4M", "4M",
         "10M", "10M", "10M", "10M", "10M", "10M","10M", "10M", "10M", "10M", "10M", "10M")
exposure <- c("BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control",
              "BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control",
              "BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control")
mouse <-c(1:6,7:12,1:6,7:12,1:6,7:12) #paired mouse covariate
sex <- c("M", "M", "M", "M", "M", "M","F", "F", "F", "F", "F", "F",
         "M", "M", "M", "M", "M", "M","F", "F", "F", "F", "F", "F",
         "M", "M", "M", "M", "M", "M","F", "F", "F", "F", "F", "F")
design_combined <- data.frame(n,age,exposure,mouse,sex)

########### DSS Package - multifactor modeling ###############
biocLite("DSS")
library(DSS)

#2/23/17 AND 5/12/17 - ONES server - Multifactorial models for examining age and exposure effects in mice while controlling for mouse and sex:
DMLfit_combined = DMLfit.multiFactor(bsseq_combined, design=design_combined, 
                                     formula=~exposure+age+mouse+sex+age:exposure)
#4/25/17 - Try design matrix w/o interaction term
DMLfit_combined_age = DMLfit.multiFactor(bsseq_combined, design=design_combined, 
                                     formula=~exposure+age+mouse+sex) 

#3/29/18 - Try design matrix w/o interaction term
DMLfit_combined_age = DMLfit.multiFactor(bsseq_combined, design=design_combined, 
                                         formula=~exposure+age+mouse+sex) 

head(DMLfit_combined_age)

#Now, we can use pull out the age or exposure effect. 
#It is important to note that the coef parameter is the index of the coefficient 
#to be tested for being 0. Because the model (as specified by formula in DMLfit.multiFactor) 
#includes the intercept, the age effect is the 2nd column and the exposure effect
#is the 3rd column in the design matrix, so we use coef=2 or coef=3, respectively.

#Note: Results are data frames with chromosome number, CpG site position, test statistics, p-values (from
#normal distribution), and FDR.
summary(DMLfit_combined)
#Effect of exposure:
DMLtest_combined.exp = DMLtest.multiFactor(DMLfit_combined, coef=2)

#Effect of age:
DMLtest_combined.age = DMLtest.multiFactor(DMLfit_combined, coef=3)
DMLtest_combined.age = DMLtest.multiFactor(DMLfit_combined_age, coef=3) #3/29/18

#Effect of age:exposure interaction (Not part of 4/25/17 model):
DMLtest_combined.age_exp = DMLtest.multiFactor(DMLfit_combined, coef=6)

##SORT BY AGE TERM
#sort the data frames (and remove NAs) using following codes:
DMLtest_combined.age_sort1 <- DMLtest_combined.age[order(DMLtest_combined.age$pvals, na.last=NA), ]
head(DMLtest_combined.age_sort1)
#Subset data frame based on FDR cutoff of 0.05
DMLtest_combined.age_sort2 <- subset(DMLtest_combined.age_sort1, fdrs < 0.05)
head(DMLtest_combined.age_sort2) ##34688 CpG sites in data frame
summary(DMLtest_combined.age_sort2)

##3/22/18 - try DMR calling on DML age term.
dmrscombined.age <- callDMR(DMLfit_combined, p.threshold=0.05)
#Works, but ONLY at a very low p-value threshold (P=0.2 or higher)


#Write CSVs for results from models sorted by p-values.
write.csv(DMLtest_combined.age_sort2, file = "DSS_DMLtest_Age_2-4-10M Samples_042017.csv") #Updated 5/16/17 with correct order of samples
write.csv(DMLtest_combined.age_sort2, file = "DSS_DMLtest_Age_2-4-10M Samples_042517.csv")
write.csv(DMLtest_combined.age_sort2, file = "DSS_DMLtest_Age_2-4-10M Samples_060717.csv") #Reran DMLtest on 6/7/17 to ensure correct variable designations
write.csv(DMLtest_combined.age_sort2, file = "DSS_DMLtest_Age_2-4-10M Samples_032918.csv") #Reran DMLtest on 6/7/17 to ensure correct variable designations


##SORT BY EXPOSURE TERM
#sort the data frames (and remove NAs) using following codes:
DMLtest_combined.exp_sort1 <- DMLtest_combined.exp[order(DMLtest_combined.exp$pvals, na.last=NA), ]
head(DMLtest_combined.exp_sort1)
#Subset data frame based on FDR cutoff of 0.05
DMLtest_combined.exp_sort2 <- subset(DMLtest_combined.exp_sort1, fdrs < 0.05)
head(DMLtest_combined.exp_sort2) ##36497 CpG sites in data frame
summary(DMLtest_combined.exp_sort2)

#Write CSVs for exposure term p-values from models.
write.csv(DMLtest_combined.exp_sort2, file = "DSS_DMLtest_Exposure_2-4-10M Samples_042017.csv") #Updated 5/16/17 with correct order of samples
write.csv(DMLtest_combined.exp_sort2, file = "DSS_DMLtest_Exposure_2-4-10M Samples_042517.csv")
write.csv(DMLtest_combined.exp_sort2, file = "DSS_DMLtest_Exposure_2-4-10M Samples_060717.csv") #Reran DMLtest on 6/7/17 to ensure correct variable designations

##SORT BY AGE:EXPOSURE INTERACTION TERM (NOT ON 4/25/17)
#sort the data frames (and remove NAs) using following codes:
DMLtest_combined.age_exp_sort1 <- DMLtest_combined.age_exp[order(DMLtest_combined.age_exp$pvals, na.last=NA), ]
head(DMLtest_combined.age_exp_sort1)
#Subset data frame based on FDR cutoff of 0.05
DMLtest_combined.age_exp_sort2 <- subset(DMLtest_combined.age_exp_sort1, fdrs < 0.05)
head(DMLtest_combined.age_exp_sort2) ##33566 CpG sites in data frame
summary(DMLtest_combined.age_exp_sort2)
#Write CSVs for male and female age:exposure interaction term p-values from models.
write.csv(DMLtest_combined.age_exp_sort2, file = "DSS_DMLtest_AgeExp Interaction_2-4-10M Samples_042017.csv") #Updated 5/16/17 with correct order of samples
write.csv(DMLtest_combined.age_exp_sort2, file = "DSS_DMLtest_AgeExp Interaction_2-4-10M Samples_060717.csv") #Updated 5/16/17 with correct order of samples


#Finally, look at the distributions of p-values
#figure showing the distributions of test statistics and p-values from this example dataset
par(mfrow=c(1,2))
hist(DMLtest_combined.age_sort1$stat, 50, main="Combined Age test statistics", xlab="")
hist(DMLtest_combined.age_sort1$pvals, 50, main="Combined Age P-values", xlab="")

hist(DMLtest_combined.exp_sort1$stat, 50, main="Combined Exp test statistics", xlab="")
hist(DMLtest_combined.age_sort1$pvals, 50, main="Combined Exp P-values", xlab="")

hist(DMLtest_combined.age_exp_sort1$stat, 50, main="Combined Age:Exp test statistics", xlab="")
hist(DMLtest_combined.age_exp_sort1$pvals, 50, main="Combined Age:Exp P-values", xlab="")

################ IMPORTANT NOTE: ###################
##To obtain a workable form of multifactorial results for DMR calling, merge multifactorial 
##CpG sites with those found using simple DMLtest function:

makeBSseqData(dat, sampleNames)


#DMLtest on combined 2M, 4M, and 10M data - ONES Server 4/20/17
DMLtest_age <- DMLtest(bsseq_combined, c("2M_1","2M_2","2M_3","2M_4","2M_5","2M_6", "2M_7", "2M_8","2M_9","2M_10","2M_11","2M_12"), 
                                           c("4M_1", "4M_2","4M_3","4M_4","4M_5","4M_6","4M_7", "4M_8","4M_9","4M_10","4M_11","4M_12"),
                                           c("10M_1","10M_2","10M_3","10M_4","10M_5","10M_6", "10M_7","10M_8","10M_9","10M_10","10M_11", "10M_12"), 
                                           equal.disp = FALSE, smoothing = FALSE)
head(DMLtest_age)

#Note: Exposure is coded as follows
#exposure <- c("BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control",
#              "BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control",
#              "BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control","BPA", "Control")

#8/2/17 -- Created new reference DMLtest spreadsheet for BPA vs. control comparison
#To be used when calling BPA-related DMRs from DMLmultifactor fit results (which must be merged)
DMLtest_exposure <- DMLtest(bsseq_combined, c("2M_1","2M_3","2M_5","2M_7","2M_9","2M_11", 
                                              "4M_1","4M_3","4M_5","4M_7","4M_9","4M_11",
                                              "10M_1","10M_3","10M_5","10M_7","10M_9","10M_11"), 
                                            c("2M_2","2M_4","2M_6","2M_8","2M_10","2M_12", 
                                              "4M_2","4M_4","4M_6","4M_8","4M_10","4M_12",
                                              "10M_2","10M_4","10M_6","10M_8","10M_10","10M_12"), 
                       equal.disp = FALSE, smoothing = FALSE)
head(DMLtest_exposure)


#Read in multifactorial results from .csv files:
DMLtest_combined.age_sort2 <- read.csv("DSS_DMLtest_Age_2-4-10M Samples_042017.csv") #Updated 5/16/17 with correct order of samples
DMLtest_combined.age_sort2 <- read.csv("DSS_DMLtest_Age_2-4-10M Samples_042517.csv") #model w/o interaction term 
DMLtest_combined.age_sort2 <- read.csv("DSS_DMLtest_Age_2-4-10M Samples_060717.csv") #6/7/17 - model w/ interaction term 

head(DMLtest_combined.age_sort2)
DMLtest_combined.exp_sort2 <- read.csv("DSS_DMLtest_Exposure_2-4-10M Samples_042017.csv") #Updated 5/16/17 with correct order of samples
DMLtest_combined.exp_sort2 <- read.csv("DSS_DMLtest_Exposure_2-4-10M Samples_042517.csv")#model w/o interaction term 
DMLtest_combined.exp_sort2 <- read.csv("DSS_DMLtest_Exposure_2-4-10M Samples_060717.csv")#6/7/17 - model w/ interaction term 

DMLtest_combined.age_exp_sort2 <- read.csv("DSS_DMLtest_AgeExp Interaction_2-4-10M Samples_042017.csv") #Updated 5/16/17 with correct order of samples
DMLtest_combined.age_exp_sort2 <- read.csv("DSS_DMLtest_AgeExp Interaction_2-4-10M Samples_060717.csv") #6/7/17 - model w/ interaction and correct order of samples

#Merge data frames based on shared chr, pos variables
#Note: this first merge will ONLY include common cases in both datasets
#Updated on 8/2/17:
DMLtest_combined.age_merged <- merge(DMLtest_age, DMLtest_combined.age_sort2, by=c("chr","pos")) 
DMLtest_combined.exp_merged <- merge(DMLtest_exposure, DMLtest_combined.exp_sort2, by=c("chr","pos")) 
DMLtest_combined.age_exp_merged <- merge(DMLtest_age, DMLtest_combined.age_exp_sort2, by=c("chr","pos")) 

##Write .csv files of merged DMLtest_combined.age results
write.csv(DMLtest_combined.age_merged, file = "DSS_DMLtest_Age_Combined_MergedData_042017.csv") #Updated 5/16/17 with correct order of samples
write.csv(DMLtest_combined.age_merged, file = "DSS_DMLtest_Age_Combined_MergedData_042517.csv")
write.csv(DMLtest_combined.age_merged, file = "DSS_DMLtest_Age_Combined_MergedData_060717.csv")

write.csv(DMLtest_combined.exp_merged, file = "DSS_DMLtest_Exposure_Combined_MergedData_042017.csv") #Updated 5/16/17 with correct order of samples
write.csv(DMLtest_combined.exp_merged, file = "DSS_DMLtest_Exposure_Combined_MergedData_042517.csv")
write.csv(DMLtest_combined.exp_merged, file = "DSS_DMLtest_Exposure_Combined_MergedData_060717.csv")
write.csv(DMLtest_combined.exp_merged, file = "DSS_DMLtest_Exposure_Combined_MergedData_080217.csv") #8/2/17 -- Updated DMLtest list to reflect BPA-related baseline comparison


write.csv(DMLtest_combined.age_exp_merged, file = "DSS_DMLtest_Age_Exp_Combined_MergedData_042017.csv") #Updated 5/16/17 with correct order of samples
write.csv(DMLtest_combined.age_exp_merged, file = "DSS_DMLtest_Age_Exp_Combined_MergedData_060717.csv") #Updated 5/16/17 with correct order of samples

###############################################################################
#8/1/17 -- Subset data frame based on difference in methylation cutoff of 0.10
#DOES NOT WORK in DSS DMLfit.multifactor because NO "diff" variable in multifactor modeling.
#However, this is available in the "merged" datasets created above!

##Load in merged datasets
setwd("C:/Users/jjkoch/Dropbox/Next-gen_Sequencing_Data/DSS DM Outputs/2-4-10M_DMLs_060717")
DMLtest_combined.age_merged <- read.csv("DSS_DMLtest_Age_Combined_MergedData_060717.csv")
DMLtest_combined.age_exp_merged <- read.csv("DSS_DMLtest_Age_Exp_Combined_MergedData_060717.csv")

setwd("C:/Users/jjkoch/Dropbox/Next-gen_Sequencing_Data/DSS DM Outputs/2-4-10M_DMLs_080217")
DMLtest_combined.exp_merged <- read.csv("DSS_DMLtest_Exposure_Combined_MergedData_080217.csv")

#Institute more stringent cutoffs for DMC-calling
#5% change in methylation cutoff for age effects
DMLtest_combined.age_merged_sort1 <- subset(DMLtest_combined.age_merged, abs(diff) > 0.05)
#28196 DMCs

#5% change in methylation cutoff for exposure effects
DMLtest_combined.exp_merged_sort1 <- subset(DMLtest_combined.exp_merged, abs(diff) > 0.05)
#50330 DMCs

#5% change in methylation cutoff for age effects
DMLtest_combined.age_exp_merged_sort1 <- subset(DMLtest_combined.age_exp_merged, abs(diff) > 0.05)
#37511 DMCs

write.csv(DMLtest_combined.age_merged_sort1, file = "DSS_Age_2-4-10M_MergedData_Stringent_DMCs_080117.csv")
write.csv(DMLtest_combined.exp_merged_sort1, file = "DSS_Exposure_2-4-10M_MergedData_Stringent_DMCs_080217.csv")
write.csv(DMLtest_combined.age_exp_merged_sort1, file = "DSS_AgeExp_2-4-10M_MergedData_Stringent_DMCs_080117.csv")

#Call DMRs on DMLtest_combined.age_merged results!

##In Excel, truncate the DMLtest_combined.age_merged data.frame to remove the last four
##columns from the multifactorial model.

#Age
DMLtest_combined.age_merged1 <- read.csv("DSS_DMLtest_Age_Combined_MergedData_Truncated_042017.csv")
DMLtest_combined.age_merged1 <- read.csv("DSS_DMLtest_Age_Combined_MergedData_Truncated_060717.csv")
#More stringent list of DMCs (based on >5% methylation change):
DMLtest_combined.age_merged1 <- read.csv("DSS_Age_2-4-10M_MergedData_Stringent_DMCs_Truncated_080117.csv")

#BPA exposure
DMLtest_combined.exp_merged1 <- read.csv("DSS_DMLtest_Exposure_Combined_MergedData_Truncated_042017.csv")
DMLtest_combined.exp_merged1 <- read.csv("DSS_DMLtest_Exposure_Combined_MergedData_Truncated_060717.csv")
DMLtest_combined.exp_merged1 <- read.csv("DSS_DMLtest_Exposure_Combined_MergedData_Truncated_080217.csv")
#More stringent list of DMCs (based on >5% methylation change):
DMLtest_combined.exp_merged1 <- read.csv("DSS_Exposure_2-4-10M_MergedData_Stringent_DMCs_Truncated_080217.csv")

#Age:BPA interaction
DMLtest_combined.age_exp_merged1 <- read.csv("DSS_DMLtest_Age_Exp_Combined_MergedData_Truncated_042017.csv")
DMLtest_combined.age_exp_merged1 <- read.csv("DSS_DMLtest_Age_Exp_Combined_MergedData_Truncated_060717.csv")
#More stringent list of DMCs (based on >5% methylation change):
DMLtest_combined.age_exp_merged1 <- read.csv("DSS_AgeExp_2-4-10M_MergedData_Stringent_DMCs_080117.csv")

#callDMR for DMR calling:
#Default settings - callDMR(DMLresult, delta=0, p.threshold=1e-5,
#minlen=50, minCG=3, dis.merge=100, pct.sig=0.5) #see DSS vignette for details on default parameters!


dmrscombined.age <- callDMR(DMLtest_combined.age_merged1, delta=0.05, p.threshold=0.05)
## 6 DMRs using p-value cutoff of 1E-5; 25 DMRs for p-value cutoff of 1E-3; 55 DMRs using p-value cutoff of 1E-2.
##152 DMRs called using p-value cutoff of 0.05
#5/16/17 -- Updated BPA/Control designation -- 141 DMRs called using p-value cutoff of 0.05
#8/2/17 -- 20 DMRs called using p-value cutoff of 0.05 and stringent DMC calling threshold (>5%) 

dmrscombined.exp <- callDMR(DMLtest_combined.exp_merged1, delta=0.05, p.threshold=0.05)
## 4 DMRs using p-value cutoff of 1E-5; 5 DMRs for p-value cutoff of 1E-3; 8 DMRs using p-value cutoff of 1E-2.
##5/16/17 -- 17 DMRs called using p-value cutoff of 0.05
##8/2/17 --21 DMRs called using updated BPA-related baseline comparison for data merging
##8/2/17 -- 22 DMRs called using updated BPA-related baseline comparison for data merging and 
#stringent (>5% methylation change) DMC-calling method.

dmrscombined.age_exp <- callDMR(DMLtest_combined.age_exp_merged1, delta=0.05, p.threshold=0.05)
## 0 DMRs using p-value cutoff of 1E-5; 1 DMRs for p-value cutoff of 1E-3; 2 DMRs using p-value cutoff of 1E-2.
##15 DMRs called using p-value cutoff of 0.05
##5/16/17 -- 1 DMR called using p-value cutoff of 0.05
#8/2/17 -- 0 DMRs called using p-value cutoff of 0.05 and stringent DMC calling threshold (>10%) 

#View top hits from each set of called DMRs
head(dmrscombined.age) 
head(dmrscombined.exp)
head(dmrscombined.age_exp)

##Write .csv of results!
write.csv(dmrscombined.age, file = "DSS_Age_2-4-10M_MergedData_Stringent_DMRs_080217.csv")
write.csv(dmrscombined.age, file = "DSS_Age_2-4-10M_MergedData_DMRs_042017.csv")
write.csv(dmrscombined.exp, file = "DSS_Exposure_2-4-10M_MergedData_DMRs_080217.csv")
write.csv(dmrscombined.exp, file = "DSS_Exposure_2-4-10M_MergedData_Stringent_DMRs_080217.csv")
write.csv(dmrscombined.age_exp, file = "DSS_AgeExp_2-4-10M_MergedData_DMRs_042017.csv")
##View a single DMR from the list:
windows() ##This creates new window for plotting
#x11() #Creates new window in Linux
showOneDMR(dmrscombined.age[1,], bsseq_combined) #DOES NOT WORK FOR COMBINED DATA
##PLOT IS TOO LARGE FOR PLOTTING FRAME.

#Conclusion:Small number of DMRs by age when adjusting for paired mouse effect, exposure, and age:exposure 
#interaction.Interestingly, largest # of DML are seen for age:exp interaction, but this
#set of DMLs also has the fewest called DMRs. Maybe there is environmental deflection at 
#single CpG sites throughout the genome that is not reflected in DMR calling method 
#(Not affecting entire regions?!)


###############BED format:
#In excel, remove top row of variable names, then save file as .bed file.
#This BED file can then be read in by annotatr to annotate gene regions

###### 2/27/17 - Annotate DMRs using annotatr #######

##Annotate CpG sites using the annotatr package
source("https://bioconductor.org/biocLite.R")
biocLite("annotatr")
library(annotatr)
biocLite("regioneR")
library(regioneR)

setwd("C:/Users/jjkoch/Dropbox/Next-gen_Sequencing_Data/DSS DM Outputs/2-4-10M_DMRs_042017")

dmrs_combined_age = read_regions("DSS_Age_2-4-10M_MergedData_DMRs_BED format_042017.bed", genome = 'mm10')
dmrs_combined_exp = read_regions("DSS_Exposure_2-4-10M_MergedData_DMRs_BED format_042017.bed", genome = 'mm10')
dmrs_combined_age.exp = read_regions("DSS_AgeExp_2-4-10M_MergedData_DMRs_BED format_042017.bed", genome = 'mm10')

print(dmrs_combined_age) #141 ranges
print(dmrs_combined_exp) #17 ranges
print(dmrs_combined_age.exp) #1 range

#8/2/17 -- Read in DMRs called using more stringent (>5%) cutoff for DMCs
setwd("C:/Users/jjkoch/Dropbox/Next-gen_Sequencing_Data/DSS DM Outputs/2-4-10M_DMRs_080217")

dmrs_combined_age = read_regions("DSS_Age_2-4-10M_MergedData_Stringent_DMRs_BEDformat_080217.bed", genome = 'mm10')
dmrs_combined_exp = read_regions("DSS_Exposure_2-4-10M_MergedData_Stringent_DMRs_BEDformat_080217.bed", genome = 'mm10')

print(dmrs_combined_age) #20 ranges
print(dmrs_combined_exp) #22 ranges
#print(dmrs_combined_age.exp) #0 ranges

# Select and build annotations for mm10
annots = c('mm10_basicgenes','mm10_cpgs')
annotations = build_annotations(genome = 'mm10', annotations = annots)

##Now annotate!
dmrs_combined_age_annotated = annotate_regions(regions = dmrs_combined_age, annotations = annotations, ignore.strand = TRUE)
print(dmrs_combined_age_annotated)
dmrs_combined_exp_annotated = annotate_regions(regions = dmrs_combined_exp, annotations = annotations, ignore.strand = TRUE)
print(dmrs_combined_exp_annotated)
dmrs_combined_age.exp_annotated = annotate_regions(regions = dmrs_combined_age.exp, annotations = annotations, ignore.strand = TRUE)
print(dmrs_combined_age.exp_annotated)

# Note: annotatr output is a special type of dataframe using the dplyr package.
# By default, dplyr::tbl_df objects have nice printing properties, but it
# hides extra columns that would ordinarily wrap. You can see them all with:
head(as.data.frame(dmrs_combined_age_annotated))
head(as.data.frame(dmrs_combined_exp_annotated))
head(as.data.frame(dmrs_combined_age.exp_annotated))
## You could write the gene region annotations to a tab-delimited text file with:
write.table(as.data.frame(dmrs_combined_age_annotated),file="Mouse_Blood_2-4-10M DMRs_Age_annotations_042017.txt",sep="\t",quote=F)
write.table(as.data.frame(dmrs_combined_age_annotated),file="Mouse_Blood_2-4-10M DMRs_Age_Stringent_annotations_080217.txt",sep="\t",quote=F)

write.table(as.data.frame(dmrs_combined_exp_annotated),file="Mouse_Blood_2-4-10M DMRs_Exposure_annotations_042017.txt",sep="\t",quote=F)
write.table(as.data.frame(dmrs_combined_exp_annotated),file="Mouse_Blood_2-4-10M DMRs_Exposure_Stringent_annotations_080217.txt",sep="\t",quote=F)

write.table(as.data.frame(dmrs_combined_age.exp_annotated),file="Mouse_Blood_2-4-10M DMRs_AgeExp_annotations_042017.txt",sep="\t",quote=F)


##It's also important to create a random annotation file for comparison
# Randomize the input regions - Updated method on 5/16/17 that only considers
# ERRBS data as input for randomization.
source("https://bioconductor.org/biocLite.R")
biocLite("regioneR")
library(regioneR) 
#Examine bsseq data imported earlier in R session:
bsseq_combined
head(granges(bsseq_combined))

bsseq_universe <- granges(bsseq_combined)
#universe -- ERRBS input data; generally, a region set in any of the formats accepted by toGRanges (GenomicRanges,
#data.frame, etc...)

dm_random_age_combined = resampleRegions(
  A = dmrs_combined_age, 
  universe = bsseq_universe,
  allow.overlaps = TRUE, 
  per.chromosome = TRUE) #Generates random distribution from ERRBS data universe

dm_random_exp_combined = resampleRegions(
  A = dmrs_combined_exp, 
  universe = bsseq_universe,
  allow.overlaps = TRUE, 
  per.chromosome = TRUE)#Generates random distribution from ERRBS data universe

dm_random_age.exp_combined = resampleRegions(
  A = dmrs_combined_age.exp, 
  universe = bsseq_universe,
  allow.overlaps = TRUE, 
  per.chromosome = TRUE)#Generates random distribution from ERRBS data universe

############# OLD METHOD ##################################################
##### Random distributions from entire genome using annotatr
#dm_random_age_combined = randomize_regions(
#  regions = dmrs_combined_age,
#  allow.overlaps = TRUE,
#  per.chromosome = TRUE) #Generates random distribution from entire genome

#dm_random_exp_combined = randomize_regions(
#  regions = dmrs_combined_exp,
#  allow.overlaps = TRUE,
#  per.chromosome = TRUE) #Generates random distribution from entire genome

#dm_random_age.exp_combined = randomize_regions(
#  regions = dmrs_combined_age.exp,
#  allow.overlaps = TRUE,
#  per.chromosome = TRUE) #Generates random distribution from entire genome
###########################################################################

# Annotate the random regions using the same annotations as above
# These will be used in later functions
dm_random_age_combined_annotated = annotate_regions(
  regions = dm_random_age_combined,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE)

dm_random_exp_combined_annotated = annotate_regions(
  regions = dm_random_exp_combined,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE)

dm_random_age.exp_combined_annotated = annotate_regions(
  regions = dm_random_age.exp_combined,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE)

# Usage of summarize functions with defaults
##Compares annotated DMRs to random annotations
summarize_annotations(dmrs_combined_age_annotated, dm_random_age_combined_annotated, quiet = FALSE)
summarize_annotations(dmrs_combined_exp_annotated, dm_random_exp_combined_annotated, quiet = FALSE)
summarize_annotations(dmrs_combined_age.exp_annotated, dm_random_age.exp_combined_annotated, quiet = FALSE)
##Larger number of annotations in gene exons in real data compared to randomized annotations. Other categories are similar.

# Can also pull out ONLY the number of regions per annotation type
dm_annsum_age = summarize_annotations(
  annotated_regions = dmrs_combined_age_annotated,
  quiet = TRUE)
print(dm_annsum_age)

dm_annsum_exp = summarize_annotations(
  annotated_regions = dmrs_combined_exp_annotated,
  quiet = TRUE)
print(dm_annsum_exp)

dm_annsum_age.exp = summarize_annotations(
  annotated_regions = dmrs_combined_age.exp_annotated,
  quiet = TRUE)
print(dm_annsum_age.exp)
##Plotting Annotation results:
# View the number of regions per annotation and include the annotation
# of randomized regions
annots_order = c(
  'mm10_cpg_inter',
  'mm10_cpg_islands',
  'mm10_cpg_shores',
  'mm10_cpg_shelves',
  'mm10_genes_promoters',
  'mm10_genes_3UTRs',
  'mm10_genes_5UTRs',
  'mm10_genes_exons',
  'mm10_genes_introns')
dmrs_age_annotations_wrandom = plot_annotation(
  annotated_regions = dmrs_combined_age_annotated,
  annotated_random = dm_random_age_combined_annotated,
  annotation_order = annots_order,
  plot_title = 'Dist. of DMRs by Age (with rndm.) ',
  x_label = 'Annotations',
  y_label = 'Count')
pdf("Mouse_2-4-10M_Combined_Age_Stringent_DMR Distribution_080217.pdf") #5/17/17: Updated Plot created for ALL ERRBS data
print(dmrs_age_annotations_wrandom) #Plot generated on 5/17/17; 8/2/17 - Updated, more stringent cutoff version 
dev.off()


#Repeat for exposure effects in ALL samples - 4/20/17
dmrs_exp_annotations_wrandom = plot_annotation(
  annotated_regions = dmrs_combined_exp_annotated,
  annotated_random = dm_random_exp_combined_annotated,
  annotation_order = annots_order,
  plot_title = 'Dist. of DMRs by BPA Exposure (with rndm.) ',
  x_label = 'Annotations',
  y_label = 'Count')
pdf("Mouse_2-4-10M_Combined_Exp_Stringent_DMR Distribution_080217.pdf") #5/17/17: Updated Plot created for ALL ERRBS data
print(dmrs_exp_annotations_wrandom) #Plot generated on 5/17/17; 8/2/17 - Updated, more stringent version
dev.off()

#Repeat for age:exp interaction in ALL samples - 2/27/17
dmrs_age.exp_annotations_wrandom = plot_annotation(
  annotated_regions = dmrs_combined_age.exp_annotated,
  annotated_random = dm_random_age.exp_combined_annotated,
  annotation_order = annots_order,
  plot_title = 'Dist. of DMRs by Age*Exp (with rndm.) ',
  x_label = 'Annotations',
  y_label = 'Count')
pdf("Mouse_2-4-10M_Combined_AgeExp_DMR Distribution_051717.pdf") #5/17/17: Updated Plot created for ALL ERRBS data
print(dmrs_age.exp_annotations_wrandom) #Plot generated on 5/17/17
dev.off()


#####################################################
###### 5/17/17 - Annotate individual DMCs using annotatr #######
###### 8/2/17 -- Updated version of DMC annotation (stringent >5% change dataset)

##Annotate CpG sites using the annotatr package
source("https://bioconductor.org/biocLite.R")
biocLite("S4Vectors")
biocLite("RSQLite")
biocLite("annotatr")
install.packages("Matrix")
library(annotatr)
biocLite("GenomicRanges")
library(GenomicRanges)
setwd("/jjkoch/03-methylation_call")
setwd("C:/Users/jjkoch/Dropbox/Next-gen_Sequencing_Data/DSS DM Outputs/2-4-10M_DMLs_080217")

#Read in DMLs (BED format -- note: the two chromosomal location columns (start and end) have the same value.
dmls_combined_age = read_regions("DSS_DMLtest_Age_2-4-10M Samples_BED format_060717.bed", genome = 'mm10')
dmls_combined_exp = read_regions("DSS_DMLtest_Exposure_2-4-10M Samples_BED format_060717.bed", genome = 'mm10')
dmls_combined_age.exp = read_regions("DSS_DMLtest_AgeExp Interaction_2-4-10M Samples_BED format_060717.bed", genome = 'mm10')

#8/2/17 - Updated DMC annotation
dmls_combined_age = read_regions("DSS_Age_2-4-10M_MergedData_Stringent_DMCs_BEDformat_080217.bed", genome = 'mm10')
dmls_combined_exp = read_regions("DSS_Exposure_2-4-10M_MergedData_Stringent_DMCs_BEDformat_080217.bed", genome = 'mm10')
dmls_combined_age.exp = read_regions("DSS_AgeExp_2-4-10M_MergedData_Stringent_DMCs_BEDformat_080217.bed", genome = 'mm10')


print(dmls_combined_age) #28196 Ranges
print(dmls_combined_exp) #50330 Ranges
print(dmls_combined_age.exp) #37511 Ranges

#Annotatr adds 1 to the second column when organizing the data, so reset the "start" and "end" 
#columns to the same value!
start(dmls_combined_age) = end(dmls_combined_age)
start(dmls_combined_exp) = end(dmls_combined_exp)
start(dmls_combined_age.exp) = end(dmls_combined_age.exp)

# Select and build annotations for mm10
annots = c('mm10_basicgenes','mm10_cpgs')
annotations = build_annotations(genome = 'mm10', annotations = annots)

##Now annotate the DMLs to the genome
dmls_combined_age_annotated = annotate_regions(regions = dmls_combined_age, annotations = annotations, ignore.strand = TRUE)
print(dmls_combined_age_annotated) #8/2/17 - 73411 annotated Ranges
dmls_combined_exp_annotated = annotate_regions(regions = dmls_combined_exp, annotations = annotations, ignore.strand = TRUE)
print(dmls_combined_exp_annotated) #8/2/17 - 132078 annotated Ranges
dmls_combined_age.exp_annotated = annotate_regions(regions = dmls_combined_age.exp, annotations = annotations, ignore.strand = TRUE)
print(dmls_combined_age.exp_annotated) #8/2/17 - 103730 annotated Ranges

# Note: annotatr output is a special type of dataframe using the dplyr package.
# By default, dplyr::tbl_df objects have nice printing properties, but it
# hides extra columns that would ordinarily wrap. You can see them all with:
head(as.data.frame(dmls_combined_age_annotated))
head(as.data.frame(dmls_combined_exp_annotated))
head(as.data.frame(dmls_combined_age.exp_annotated))
## You could write the gene region annotations to a tab-delimited text file with:
write.table(as.data.frame(dmls_combined_age_annotated),file="Mouse_Blood_2-4-10M DMLs_Age_annotations_051717.txt",sep="\t",quote=F)
write.table(as.data.frame(dmls_combined_exp_annotated),file="Mouse_Blood_2-4-10M DMLs_Exposure_annotations_051717.txt",sep="\t",quote=F)
write.table(as.data.frame(dmls_combined_age.exp_annotated),file="Mouse_Blood_2-4-10M DMLs_AgeExp_annotations_051717.txt",sep="\t",quote=F)

write.table(as.data.frame(dmls_combined_age_annotated),file="Mouse_Blood_2-4-10M DMLs_Age_Stringent_annotations_080217.txt",sep="\t",quote=F)
write.table(as.data.frame(dmls_combined_exp_annotated),file="Mouse_Blood_2-4-10M DMLs_Exposure_Stringent_annotations_080217.txt",sep="\t",quote=F)
write.table(as.data.frame(dmls_combined_age.exp_annotated),file="Mouse_Blood_2-4-10M DMLs_AgeExp_Stringent_annotations_080217.txt",sep="\t",quote=F)

#### 3/22/18 ####
####Use plot categorical to plot annotations
dm_order = c(
  'hyper',
  'hypo')

cpg_order = c(
  'hg19_cpg_islands',
  'hg19_cpg_shores',
  'hg19_cpg_shelves',
  'hg19_cpg_inter')

dm_vn = plot_categorical(
  annotated_regions = dmls_combined_age_annotated,
  x = 'DM_status',
  fill = 'annot.type',
  x_order = dm_order,
  fill_order = cpg_order,
  position = 'fill',
  legend_title = 'knownGene Annotations',
  x_label = 'DM status',
  y_label = 'Proportion')

##It's also important to create a random annotation file for comparison
# Randomize the input regions
dml_random_age_combined = resampleRegions(
  A = dmls_combined_age, 
  universe = bsseq_universe,
  allow.overlaps = TRUE, 
  per.chromosome = TRUE) #Generates random distribution from ERRBS data universe

dml_random_exp_combined = resampleRegions(
  A = dmls_combined_exp, 
  universe = bsseq_universe,
  allow.overlaps = TRUE, 
  per.chromosome = TRUE)#Generates random distribution from ERRBS data universe

dml_random_age.exp_combined = resampleRegions(
  A = dmls_combined_age.exp, 
  universe = bsseq_universe,
  allow.overlaps = TRUE, 
  per.chromosome = TRUE)#Generates random distribution from ERRBS data universe

# Annotate the random regions using the same annotations as above
# These will be used in later functions
dml_random_age_combined_annotated = annotate_regions(
  regions = dml_random_age_combined,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE)

dml_random_exp_combined_annotated = annotate_regions(
  regions = dml_random_exp_combined,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE)

dml_random_age.exp_combined_annotated = annotate_regions(
  regions = dml_random_age.exp_combined,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE)

# Usage of summarize functions with defaults
##Compares annotated DMLs to random annotations
summarize_annotations(dmls_combined_age_annotated, dml_random_age_combined_annotated, quiet = FALSE)
summarize_annotations(dmls_combined_exp_annotated, dml_random_exp_combined_annotated, quiet = FALSE)
summarize_annotations(dmls_combined_age.exp_annotated, dml_random_age.exp_combined_annotated, quiet = FALSE)
##Larger number of annotations in gene exons in real data compared to randomized annotations. Other categories are similar.

# Can also pull out ONLY the number of regions per annotation type
dml_annsum_age = summarize_annotations(
  annotated_regions = dmls_combined_age_annotated,
  quiet = TRUE)
print(dml_annsum_age)

dml_annsum_exp = summarize_annotations(
  annotated_regions = dmls_combined_exp_annotated,
  quiet = TRUE)
print(dml_annsum_exp)

dml_annsum_age.exp = summarize_annotations(
  annotated_regions = dmls_combined_age.exp_annotated,
  quiet = TRUE)
print(dml_annsum_age.exp)
##Plotting Annotation results:
# View the number of regions per annotation and include the annotation
# of randomized regions
annots_order = c(
  'mm10_cpg_inter',
  'mm10_cpg_islands',
  'mm10_cpg_shores',
  'mm10_cpg_shelves',
  'mm10_genes_promoters',
  'mm10_genes_3UTRs',
  'mm10_genes_5UTRs',
  'mm10_genes_exons',
  'mm10_genes_introns')
dmls_age_annotations_wrandom = plot_annotation(
  annotated_regions = dmls_combined_age_annotated,
  annotated_random = dml_random_age_combined_annotated,
  annotation_order = annots_order,
  plot_title = 'Dist. of DMCs by Age (with rndm.) ',
  x_label = 'Annotations',
  y_label = 'Count')
pdf("Mouse_2-4-10M_Combined_Age_Stringent_DML Distribution_080217.pdf") #5/17/17: Updated Plot created for ALL ERRBS data
print(dmls_age_annotations_wrandom) #Plot generated on 5/17/17; updated on 8/2/17
dev.off()


#Repeat for exposure effects in ALL samples - 5/17/17
dmls_exp_annotations_wrandom = plot_annotation(
  annotated_regions = dmls_combined_exp_annotated,
  annotated_random = dml_random_exp_combined_annotated,
  annotation_order = annots_order,
  plot_title = 'Dist. of DMCs by BPA Exposure (with rndm.) ',
  x_label = 'Annotations',
  y_label = 'Count')
pdf("Mouse_2-4-10M_Combined_Exp_Stringent_DML Distribution_0080217.pdf") #5/17/17: Updated Plot created for ALL ERRBS data
print(dmls_exp_annotations_wrandom) #Plot generated on 5/17/17; updated on 8/2/17
dev.off()

#Repeat for age:exp interaction in ALL samples - 5/17/17
dmls_age.exp_annotations_wrandom = plot_annotation(
  annotated_regions = dmls_combined_age.exp_annotated,
  annotated_random = dml_random_age.exp_combined_annotated,
  annotation_order = annots_order,
  plot_title = 'Dist. of DMCs by Age*Exp (with rndm.) ',
  x_label = 'Annotations',
  y_label = 'Count')
pdf("Mouse_2-4-10M_Combined_AgeExp_Stringent_DML Distribution_080217.pdf") #5/17/17: Updated Plot created for ALL ERRBS data
print(dmls_age.exp_annotations_wrandom) #Plot generated on 5/17/17; updated on 8/2/17
dev.off()


