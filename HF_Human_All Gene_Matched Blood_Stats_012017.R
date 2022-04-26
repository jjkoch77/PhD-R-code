
##Healthy Families 
##Integration of Physical activity into longitudinal models 

setwd("/Users/jjkoch/Dropbox/Healthy Families/Healthy Families Statistics") ##Set working directory to location of file of interest

BloodMeth1 <- read.csv("HF_Human_All Gene_Matched Samples_LME Formatted_100316.csv")

attach(BloodMeth1) ##Set created dataframe as working data set for models
#detach(BloodMeth1) ##Remove working dataframe

summary(BloodMeth1)

##Frequency count (N) for each mean methylation variable modeled later in this code:
install.packages("psych")
library(psych)
describeBy(LINE1_Mean_Methylation, Matched_Age)
describeBy(IGF2_Mean_Methylation, Matched_Age)
describeBy(H19_Mean_Methylation, Matched_Age)
describeBy(PPARA_Mean_Methylation, Matched_Age)
describeBy(LEP_Mean_Methylation, Matched_Age)
describeBy(ESR1_Mean_Methylation, Matched_Age)
describeBy(SREBF11_Mean_Methylation, Matched_Age)

##Descriptive statistics for each gene within individual age groups
##Split by bloodspot/LL blood categories
describeBy(LINE1_Mean_Methylation, Age_category)
describeBy(IGF2_Mean_Methylation, Age_category)
describeBy(H19_Mean_Methylation, Age_category)
describeBy(PPARA_Mean_Methylation, Age_category)
describeBy(LEP_Mean_Methylation, Age_category)
describeBy(ESR1_Site_1, Age_category)
describeBy(ESR1_Site_4, Age_category)
describeBy(SREBF1_Mean_Methylation, Age_category)

##Histograms to examine normality
par(mfrow=c(1,1)) ##Makes 2x1 frame histogram. 
hist(MVPA_perday,nclass=20) ##nclass=20 indicates that you want 20 bars. 
##Looks relatively normal
hist(MVPA_perday[Unmatched_Age==2],nclass=20) #Not normal, near uniform
hist(MVPA_perday[Unmatched_Age==3],nclass=20) #Kinda normal

hist(HEI_score[Unmatched_Age==1]) ##Relatively normal
hist(HEI_score[Unmatched_Age==2]) ##Relatively normal
hist(HEI_score[Unmatched_Age==3]) ##Relatively normal

par(mfrow=c(2,2))
hist(bmi,nclass=20) #non-normal, right skew
hist(bmipct,nclass=20) #non-normal, left skew (near uniform)
hist(bmiz,nclass=20) #non-normal, slightleft skew
hist(mombmi,nclass=20) #non-normal, right skew

hist(log(bmi),nclass=20) #more normal, but still slightly skewed right
hist(log(bmipct),nclass=20) #non-normal, left skew
hist(log(bmiz),nclass=20) #non-normal, left skew (one major negative outlier!)
hist(log(mombmi),nclass=20) #more normal, but still right skew

par(mfrow=c(1,1))
#hist(log(Mean.Methylation),nclass=20) #Checking whether log-transformation helps with normality. In this case, it looks fairly normal before transformation.
##Note: Repeated measures, so log transormation is probably not appropriate here.

##Scatterplots to compare BMIz to MVPA_perday
plot(bmiz,MVPA_perday) ##Scatterplot phenotype vs. MVPA per day 
plot(bmiz[Unmatched_Age==2],MVPA_perday[Unmatched_Age==2]) ##Scatterplot phenotype vs. activity levels for age group 2
plot(bmiz[Unmatched_Age==3],MVPA_perday[Unmatched_Age==3]) ##Scatterplot phenotype vs. activity levels for age group 3

##Scatterplots to compare 
plot(HEI_score,bmiz) ##Scatterplot phenotype vs. MVPA per day 
plot(HEI_score[Unmatched_Age==1],bmiz[Unmatched_Age==1]) ##Scatterplot phenotype vs. dietary recall for age group 3
abline(lm(bmiz[Unmatched_Age==1] ~ HEI_score[Unmatched_Age==1])) #BMIz decreases with increasing HEI_score
plot(HEI_score[Unmatched_Age==2],bmiz[Unmatched_Age==2]) ##Scatterplot phenotype vs. dietary recall for age group 3
abline(lm(bmiz[Unmatched_Age==2] ~ HEI_score[Unmatched_Age==2])) #BMIz decreases with increasing HEI_score
plot(HEI_score[Unmatched_Age==3],bmiz[Unmatched_Age==3]) ##Scatterplot phenotype vs. dietary recall for age group 3
abline(lm(bmiz[Unmatched_Age==3] ~ HEI_score[Unmatched_Age==3])) #BMIz slightly increases with increasing HEI_score

##Check statistical significance of relationships between HEI_score and BMIz
HEI_model0<-lm(bmiz ~ HEI_score + Unmatched_Age + childgender)
summary(HEI_model0) # p>0.232

##Create separate versions of data for each age group
BloodMeth1_age1 <- BloodMeth1[which(BloodMeth1$Unmatched_Age==1),]
BloodMeth1_age2 <- BloodMeth1[which(BloodMeth1$Unmatched_Age==2),]
BloodMeth1_age3 <- BloodMeth1[which(BloodMeth1$Unmatched_Age==3),]

##Separate models for each age group:
HEI_model1<-lm(bmiz ~ HEI_score, data=BloodMeth1_age1)
summary(HEI_model1) ##Beta = -0.05; p=0.0545 -- Warrants further examination!
HEI_model2<-lm(bmiz ~ HEI_score, data=BloodMeth1_age2)
summary(HEI_model2) ##Beta = -0.022; p=0.167
HEI_model3<-lm(bmiz ~ HEI_score, data=BloodMeth1_age3)
summary(HEI_model3) ##Beta = 0.03; p=0.399

##Add gender as covariate
HEI_model1a<-lm(bmiz ~ HEI_score + childgender, data=BloodMeth1_age1)
summary(HEI_model1a) ## HEI Beta = -0.046 and p=0.0768; gender p=0.323 
HEI_model2a<-lm(bmiz ~ HEI_score + childgender, data=BloodMeth1_age2)
summary(HEI_model2a) ##HEI p=0.599; gender p=0.00309
HEI_model3a<-lm(bmiz ~ HEI_score + childgender, data=BloodMeth1_age3)
summary(HEI_model3a) ##HEI p=0.207; gender p=0.064

##Add MVPA_perday as covariate
HEI_model1b<-lm(bmiz ~ HEI_score + childgender + MVPA_perday, data=BloodMeth1_age1)
summary(HEI_model1b) ##No MVPA_perday for group 1; no results!
HEI_model2b<-lm(bmiz ~ HEI_score + childgender + MVPA_perday, data=BloodMeth1_age2)
summary(HEI_model2b) ##HEI p=0.663; gender p=0.00634; MVPA_perday p=0.27228
HEI_model3b<-lm(bmiz ~ HEI_score + childgender + MVPA_perday, data=BloodMeth1_age3)
summary(HEI_model3b) ##HEI p=0.290; gender p=0.269; MVPA_perday p=0.282

#Split by sexes (1 vs. 2)
boxplot(MVPA_perday ~ childgender, main = "Absolute MVPA per day by Sex", 
        ylab = "MVPA (counts/day)", las = 1)
##Seems to be a subtle difference of ~10-20 avg MVPA counts per day between sexes.
boxplot(HEI_score ~ childgender, main = "HEI score by Sex", 
        ylab = "Average HEI score", las = 1)
##Seems to be small difference of ~2-3 HEI by sex.

#Split Age group 2 by sexes (1 vs. 2)
boxplot(MVPA_perday[Unmatched_Age==2] ~ childgender[Unmatched_Age==2], main = "Age Group 2 MVPA per day by Sex", 
        ylab = "MVPA (counts/day)", las = 1)
##Seems to be a subtle difference of ~5 avg MVPA counts per day between sexes.
boxplot(HEI_score[Unmatched_Age==2] ~ childgender[Unmatched_Age==2], main = "Age Group 2 HEI by Sex", 
        ylab = "Average HEI score", las = 1)
##Larger difference (~8-10 HEI) in age group 2 HEI by sex

#Split Age group 3 by sexes (1 vs. 2)
boxplot(MVPA_perday[Unmatched_Age==3] ~ childgender[Unmatched_Age==3], main = "Age Group 3 MVPA per day by Sex", 
        ylab = "MVPA (counts/day)", las = 1)
##Seems to be a larger difference of ~15 avg MVPA counts per day between sexes.
boxplot(HEI_score[Unmatched_Age==3] ~ childgender[Unmatched_Age==3], main = "Age Group 3 HEI by Sex", 
        ylab = "Average HEI score", las = 1)
##Medium-sized difference (~5 HEI) in age group 3 HEI by sex

##let's looks at the means across the different groupings:
tapply(MVPA_perday, childgender, mean, na.rm=TRUE) 
tapply(HEI_score, childgender, mean, na.rm=TRUE) 
tapply(HEI_score[Unmatched_Age==1], childgender[Unmatched_Age==1], mean, na.rm=TRUE) #shows means of methylation for different ethnicities

tapply(MVPA_perday[Unmatched_Age==2], childgender[Unmatched_Age==2], mean, na.rm=TRUE) #shows means of methylation for different ethnicities
tapply(HEI_score[Unmatched_Age==2], childgender[Unmatched_Age==2], mean, na.rm=TRUE) #shows means of methylation for different ethnicities

tapply(MVPA_perday[Unmatched_Age==3], childgender[Unmatched_Age==3], mean, na.rm=TRUE) #shows means of methylation for different races
tapply(HEI_score[Unmatched_Age==3], childgender[Unmatched_Age==3], mean, na.rm=TRUE) #shows means of methylation for different ethnicities



######################LME Models - Crude##############################

## Let's do the analysis using lmer from lme4 
install.packages("lme4")
install.packages("lmerTest")
library(lme4) 
library(lmerTest) ##Loading this package adds p-values to summary of lme4 output!! DO THIS.

##Longitudinal Modeling to test effects of age on methylation (plus covariate drift deflection).

##Given that mEsr1 showed largest effects of diet and exercise, model human ESR1 first.

##Seperate LME models for each age group.

###############Separate Longitudinal Models for Matched Age Groups###################

##Create separate databases for each Matched Age Group
AgeGroup1 <- BloodMeth1[BloodMeth1$Matched_Age==1,]

AgeGroup2 <- BloodMeth1[BloodMeth1$Matched_Age==2,]

AgeGroup3 <- BloodMeth1[BloodMeth1$Matched_Age==3,]

##Longitudinal Modeling to test effects of MVPA on mean methylation in separate matched age groups (plus covariate drift deflection).
##NOTE: MVPA data only available for two older age groups, HEI score available for ALL groups

############### LINE-1 Mean Methylation Models #############
##Model 1 - 2/2/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA available
AgeGroup1.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup1.lmer) ##HEI_score p=0.585; Age p=0.0113; BMIz p=0.4994; childgender p=0.0984

##Model 2 - 2/2/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup2.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup2.lmer) ##HEI_score p=0.00865; MVPA p-value = 0.952; Age p-value = 3.6E-06; bmiz p-value = 0.508; gender p-value = 0.103

##Model 3 - 2/2/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup3.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup3.lmer) ##HEI_score p=0.205; MVPA p-value = 0.649; Age p-value = 0.018; bmiz p-value = 0.427; gender p-value = 0.452


############### IGF2 Mean Methylation Models #############
##Model 4 - 2/2/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup4.lmer<-lmer(IGF2_Mean_Methylation~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup4.lmer) ##No terms were significant (p>0.283)

##Model 5 - 2/2/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup5.lmer<-lmer(IGF2_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup5.lmer) ##HEI p=0.198; MVPA p-value = 0.280; Age p-value = 0.000621; bmiz p-value = 0.845; gender p-value = 0.870

##Model 6 - 2/2/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup6.lmer<-lmer(IGF2_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup6.lmer) ##HEI p=0.987; MVPA p-value = 0.193; Age p-value = 0.0129; bmiz p-value = 0.197; gender p-value = 0.70851


############### H19 Mean Methylation Models #############
##Model 7 - 2/2/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup7.lmer<-lmer(H19_Mean_Methylation~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup7.lmer) ##age p-value = 0.000483; nothing else significant (p>0.508).

##Model 8 - 2/2/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup8.lmer<-lmer(H19_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup8.lmer) ##HEI p=0.116; MVPA p-value = 0.498; Age p-value = 2.22E-05; bmiz p-value = 0.406; gender p-value = 0.386

##Model 9 - 2/2/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup9.lmer<-lmer(H19_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup9.lmer) ##HEI p=0.699; MVPA p-value = 0.204; Age p-value = 0.00238; bmiz p-value = 0.614; gender p-value = 0.518


############### PPARA Mean Methylation Models #############
##Model 10 - 2/2/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup10.lmer<-lmer(PPARA_Mean_Methylation~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup10.lmer) ##age p-value = 0.00133; nothing else significant (p>0.215).

##Model 11 - 2/2/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup11.lmer<-lmer(PPARA_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup11.lmer) ##HEI p=0.925; MVPA p-value = 0.114; Age p-value = 3.95E-10; bmiz p-value = 0.356; gender p-value = 0.178

##Model 12 - 2/2/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup12.lmer<-lmer(PPARA_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup12.lmer) ##HEI p=0.509; MVPA p-value = 0.425; Age p-value = 0.00188; bmiz p-value = 0.958; gender p-value = 0.201


############### LEP Mean Methylation Models #############
##Model 13 - 2/2/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup13.lmer<-lmer(LEP_Mean_Methylation~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup13.lmer) ##HEI score p-value = 0.0927; Age p-value = 7.81E-06; bmiz p=0.8712; childgender p-value = 0.922

##Model 14 - 2/2/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup14.lmer<-lmer(LEP_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup14.lmer) ##HEI p=0.910; MVPA p-value = 0.183; Age p-value = 1.57E-07; bmiz p-value = 0.623; gender p-value = 0.929

##Model 15 - 2/2/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup15.lmer<-lmer(LEP_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup15.lmer) ##HEI p=0.286; MVPA p-value = 0.504; Age p-value = 0.0469; bmiz p-value = 0.3423; gender p-value = 0.4022


############### ESR1 Mean Methylation Models #############
##Model 16 - 2/2/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##Age p-value = 0.0497; nothing else was significant (p>0.247)

##Model 17 - 2/2/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup17.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17.lmer) ##Age p-value = 0.00391; nothing else was significant (p>0.415)

##Model 18 - 2/2/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup18.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18.lmer) ##nothing was significant (p>0.280)

##Break ESR1 down by CpG site due to differences in longitudinal patterns by site.

##Model 16a - 2/2/17 -- ESR1 Site 1 Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Site_1~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##nothing was significant (p>0.350)

##Model 16a - 2/2/17 -- ESR1 Site 4 Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Site_4~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##Age p-value = 0.00136; BMIz p-value = 0.05477

##Model 17a - 2/2/17-- ESR1 Site 1 Age Group 2; bmi z-score available for age group 2
AgeGroup17a.lmer<-lmer(ESR1_Site_1~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17a.lmer) ##nothing was significant (p>0.220)

##Model 17a - 2/2/17-- ESR1 Site 4 Age Group 2; bmi z-score available for age group 2
AgeGroup17b.lmer<-lmer(ESR1_Site_4~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17b.lmer) ##Age p-value = 0.000686; nothing else was significant (p>0.221)

##Model 18a - 2/2/17-- ESR1 Site 1 Age Group 3; bmi z-score available for age group 3
AgeGroup18a.lmer<-lmer(ESR1_Site_1~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18a.lmer) ##Age p-value = 1.08E-05; nothing else was significant (p>0.244)

##Model 18b - 2/2/17-- ESR1 Site 4 Age Group 3; bmi z-score available for age group 3
AgeGroup18b.lmer<-lmer(ESR1_Site_4~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18b.lmer) ##Age p-value = 0.0394; nothing else was significant (p>0.451)


############### SREBF1 Mean Methylation Models #############
##Model 19 - 2/2/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup19.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup19.lmer) ##HEI score p-value = 0.02053; Age p-value = 0.00497; bmiz p=0.325, gender p=0.345

##Model 20 - 2/2/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup20.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup20.lmer) ##HEI p-value = 0.546; MVPA p-value = 0.903; Age p-value = 0.000342; bmiz p-value = 0.452; gender p-value = 0.955

##Model 21 - 2/2/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup21.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup21.lmer) ##HEI p-value = 0.839; MVPA p-value = 0.921; Age p-value = 0.000124; bmiz p-value = 0.777; gender p-value = 0.177

#Repeat analyses for only SREBF1_Site 4
##Model 22 - 9/18/18 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup22.lmer<-lmer(SREBF1_Site_4~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup22.lmer) ##HEI score p-value = 0.11221; Age p-value = 0.00927; bmiz p=0.36698, gender p=0.616

##Model 23 - 9/18/18-- Age Group 2; bmi z-score available for age group 2
AgeGroup23.lmer<-lmer(SREBF1_Site_4~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup23.lmer) ##HEI p-value = 0.081491; MVPA p-value = 0.777628; Age p-value = 1.86E-05; bmiz p-value = 0.394; gender p-value = 0.943

##Model 24 - 9/18/18-- Age Group 3; bmi z-score available for age group 3
AgeGroup24.lmer<-lmer(SREBF1_Site_4~HEI_score+MVPA_perday+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup24.lmer) ##HEI p-value = 0.50561; MVPA p-value = 0.34355; Age p-value = 0.00521; bmiz p-value = 0.54977; gender p-value = 0.97028

################## Environmental Deflection #################
########################## 02/28/17 #########################

##Test for interaction between age and obesity, HEI, or MVPA for each gene##

##Linear mixed effects models separated by gene and age group:

### HEI:Age Interaction term models

############### LINE-1 Mean Methylation Models #############
##Model 1 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA available
AgeGroup1.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup1.lmer) 

##Model 2 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup2.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup2.lmer) #MVPA:Age p=0.73903

##Model 3 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup3.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup3.lmer) #MVPA:Age p=0.905


############### IGF2 Mean Methylation Models #############
##Model 4 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup4.lmer<-lmer(IGF2_Mean_Methylation~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup4.lmer) 

##Model 5 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup5.lmer<-lmer(IGF2_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup5.lmer) ##MVPA:Age p=0.782

##Model 6 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup6.lmer<-lmer(IGF2_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup6.lmer) ##MVPA:Age p=0.0934

############### H19 Mean Methylation Models #############
##Model 7 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup7.lmer<-lmer(H19_Mean_Methylation~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup7.lmer) 

##Model 8 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup8.lmer<-lmer(H19_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup8.lmer) ##MVPA:Age p=0.565

##Model 9 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup9.lmer<-lmer(H19_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup9.lmer) ##MVPA:Age p=0.91925

############### PPARA Mean Methylation Models #############
##Model 10 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup10.lmer<-lmer(PPARA_Mean_Methylation~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup10.lmer)

##Model 11 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup11.lmer<-lmer(PPARA_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup11.lmer) ##MVPA:Age p=0.6057

##Model 12 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup12.lmer<-lmer(PPARA_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup12.lmer) ##MVPA:Age p=0.510


############### LEP Mean Methylation Models #############
##Model 13 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup13.lmer<-lmer(LEP_Mean_Methylation~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup13.lmer) 

##Model 14 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup14.lmer<-lmer(LEP_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup14.lmer) ##MVPA:Age p=0.01619

##Model 15 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup15.lmer<-lmer(LEP_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup15.lmer) ##MVPA:Age p=0.0.001108


############### ESR1 Mean Methylation Models #############
##Model 16 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) 

##Model 17 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup17.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17.lmer) ##MVPA:Age p=0.480

##Model 18 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup18.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18.lmer) ##MVPA:Age p=0.637

##Break ESR1 down by CpG site due to differences in longitudinal patterns by site.

##Model 16a - 2/28/17 -- ESR1 Site 1 Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Site_1~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) 

##Model 16a - 2/28/17 -- ESR1 Site 4 Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Site_4~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) 

##Model 17a - 2/28/17-- ESR1 Site 1 Age Group 2; bmi z-score available for age group 2
AgeGroup17a.lmer<-lmer(ESR1_Site_1~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17a.lmer) ##MVPA:Age p=0.0791

##Model 17a - 2/28/17-- ESR1 Site 4 Age Group 2; bmi z-score available for age group 2
AgeGroup17b.lmer<-lmer(ESR1_Site_4~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17b.lmer) ##MVPA:Age p=0.31362

##Model 18a - 2/28/17-- ESR1 Site 1 Age Group 3; bmi z-score available for age group 3
AgeGroup18a.lmer<-lmer(ESR1_Site_1~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18a.lmer) ##MVPA:Age p=0.6952

##Model 18b - 2/28/17-- ESR1 Site 4 Age Group 3; bmi z-score available for age group 3
AgeGroup18b.lmer<-lmer(ESR1_Site_4~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18b.lmer) ##MVPA:Age p=0.239


############### SREBF1 Mean Methylation Models #############
##Model 19 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup19.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup19.lmer) 

##Model 20 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup20.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup20.lmer) ##MVPA:Age p=0.89643

##Model 21 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup21.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup21.lmer) ##MVPA:Age p=0.57614

##Model 22 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup22.lmer<-lmer(SREBF1_Site_4~HEI_score+Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup22.lmer)

##Model 23 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup23.lmer<-lmer(SREBF1_Site_4~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup23.lmer) ##Age:BMIz p=0.60420

##Model 24 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup24.lmer<-lmer(SREBF1_Site_4~HEI_score+Unmatched_Age*MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup24.lmer) ##Age:BMIz p=0.320

### HEI:Age Interaction term models ###

############### LINE-1 Mean Methylation Models #############
##Model 1 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA available
AgeGroup1.lmer<-lmer(LINE1_Mean_Methylation~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup1.lmer) ##HEI_score p=0.9005; Age p=0.391; BMIz p=0.4974; childgender p=0.0925; Age:HEI p=0.6714

##Model 2 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup2.lmer<-lmer(LINE1_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup2.lmer) ##HEI_score p=0.00331; MVPA p-value = 0.98506; Age p-value = 0.01159; bmiz p-value = 0.524; gender p-value = 0.08589
#Age:HEI p=0.12329

##Model 3 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup3.lmer<-lmer(LINE1_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup3.lmer) ##HEI_score p=0.117; MVPA p-value = 0.642; Age p-value = 0.176; bmiz p-value = 0.389; gender p-value = 0.382
#Age:HEI p =0.331


############### IGF2 Mean Methylation Models #############
##Model 4 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup4.lmer<-lmer(IGF2_Mean_Methylation~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup4.lmer) ##No terms were significant (p>0.284); Age:HEI p=0.335 

##Model 5 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup5.lmer<-lmer(IGF2_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup5.lmer) ##Age:HEI p=0.854

##Model 6 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup6.lmer<-lmer(IGF2_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup6.lmer) ##Age:HEI p=0.20927

############### H19 Mean Methylation Models #############
##Model 7 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup7.lmer<-lmer(H19_Mean_Methylation~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup7.lmer) ##Age:HEI p=0.786

##Model 8 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup8.lmer<-lmer(H19_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup8.lmer) ##Age:HEI p=0.3524

##Model 9 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup9.lmer<-lmer(H19_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup9.lmer) ##Age:HEI p=0.91925

############### PPARA Mean Methylation Models #############
##Model 10 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup10.lmer<-lmer(PPARA_Mean_Methylation~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup10.lmer) ##Age:HEI p=0.69360

##Model 11 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup11.lmer<-lmer(PPARA_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup11.lmer) ##Age:HEI p=0.42342

##Model 12 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup12.lmer<-lmer(PPARA_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup12.lmer) ##Age:HEI p=0.510


############### LEP Mean Methylation Models #############
##Model 13 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup13.lmer<-lmer(LEP_Mean_Methylation~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup13.lmer) ##Age:HEI p=0.0944

##Model 14 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup14.lmer<-lmer(LEP_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup14.lmer) ##Age:HEI p=0.40173

##Model 15 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup15.lmer<-lmer(LEP_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup15.lmer) ##Age:HEI p=0.510


############### ESR1 Mean Methylation Models #############
##Model 16 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Mean_Methylation~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##Age:HEI p=0.5958

##Model 17 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup17.lmer<-lmer(ESR1_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17.lmer) ##Age:HEI p=0.5481

##Model 18 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup18.lmer<-lmer(ESR1_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18.lmer) ##Age:HEI p=0.637

##Break ESR1 down by CpG site due to differences in longitudinal patterns by site.

##Model 16a - 2/28/17 -- ESR1 Site 1 Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Site_1~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##Age:HEI p=0.158449

##Model 16a - 2/28/17 -- ESR1 Site 4 Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Site_4~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##Age:HEI p=0.4432

##Model 17a - 2/28/17-- ESR1 Site 1 Age Group 2; bmi z-score available for age group 2
AgeGroup17a.lmer<-lmer(ESR1_Site_1~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17a.lmer) ##Age:HEI p=0.745051

##Model 17a - 2/28/17-- ESR1 Site 4 Age Group 2; bmi z-score available for age group 2
AgeGroup17b.lmer<-lmer(ESR1_Site_4~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17b.lmer) ##Age:HEI p=0.5102

##Model 18a - 2/28/17-- ESR1 Site 1 Age Group 3; bmi z-score available for age group 3
AgeGroup18a.lmer<-lmer(ESR1_Site_1~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18a.lmer) ##Age:HEI p=0.2369

##Model 18b - 2/28/17-- ESR1 Site 4 Age Group 3; bmi z-score available for age group 3
AgeGroup18b.lmer<-lmer(ESR1_Site_4~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18b.lmer) ##Age:HEI p=0.745051


############### SREBF1 Mean Methylation Models #############
##Model 19 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup19.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup19.lmer) ##Age:HEI p=0.7926

##Model 20 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup20.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup20.lmer) ##Age:HEI p=0.73815

##Model 21 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup21.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score*Unmatched_Age+MVPA_perday+bmiz+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup21.lmer) ##Age:HEI p=0.7680

##Model 22 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup22.lmer<-lmer(SREBF1_Site_4~HEI_score*Unmatched_Age+bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup22.lmer) ##Age:BMIz p=0.291

##Model 23 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup23.lmer<-lmer(SREBF1_Site_4~HEI_score*Unmatched_Age+bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup23.lmer) ##Age:BMIz p=0.667299

##Model 24 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup24.lmer<-lmer(SREBF1_Site_4~HEI_score*Unmatched_Age+bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup24.lmer) ##Age:BMIz p=0.959


######### BMIz:Age Interaction term models #########

############### LINE-1 Mean Methylation Models #############
##Model 1 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA available
AgeGroup1.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup1.lmer) ##Age:BMIz p=0.5506

##Model 2 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup2.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup2.lmer) ##Age:BMIz p=0.45627

##Model 3 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup3.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup3.lmer) ##Age:BMIz p=0.8634


############### IGF2 Mean Methylation Models #############
##Model 4 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup4.lmer<-lmer(IGF2_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup4.lmer) ##Age:BMIz p=0.708

##Model 5 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup5.lmer<-lmer(IGF2_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup5.lmer) ##Age:BMIz p=0.527897

##Model 6 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup6.lmer<-lmer(IGF2_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup6.lmer) ##Age:BMIz p=0.9496

############### H19 Mean Methylation Models #############
##Model 7 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup7.lmer<-lmer(H19_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup7.lmer) ##Age:BMIz p=0.726026

##Model 8 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup8.lmer<-lmer(H19_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup8.lmer) ##Age:BMIz p=0.806704

##Model 9 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup9.lmer<-lmer(H19_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup9.lmer) ##Age:BMIz p=0.97723

############### PPARA Mean Methylation Models #############
##Model 10 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup10.lmer<-lmer(PPARA_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup10.lmer) ##Age:BMIz p=0.59090

##Model 11 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup11.lmer<-lmer(PPARA_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup11.lmer) ##Age:BMIz p=0.42342

##Model 12 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup12.lmer<-lmer(PPARA_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup12.lmer) ##Age:BMIz p=0.510


############### LEP Mean Methylation Models #############
##Model 13 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup13.lmer<-lmer(LEP_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup13.lmer) ##Age:BMIz p=0.353

##Model 14 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup14.lmer<-lmer(LEP_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup14.lmer) ##Age:BMIz p=0.63938

##Model 15 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup15.lmer<-lmer(LEP_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup15.lmer) ##Age:BMIz p=0.7331


############### ESR1 Mean Methylation Models #############
##Model 16 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##Age:BMIz p=0.6831

##Model 17 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup17.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17.lmer) ##Age:BMIz p=0.96999

##Model 18 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup18.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18.lmer) ##Age:BMIz p=0.4699

##Break ESR1 down by CpG site due to differences in longitudinal patterns by site.

##Model 16a - 2/28/17 -- ESR1 Site 1 Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Site_1~HEI_score+Unmatched_Age*bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##Age:BMIz p=0.34265

##Model 16a - 2/28/17 -- ESR1 Site 4 Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Site_4~HEI_score+Unmatched_Age*bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##Age:BMIz p=0.4432

##Model 17a - 2/28/17-- ESR1 Site 1 Age Group 2; bmi z-score available for age group 2
AgeGroup17a.lmer<-lmer(ESR1_Site_1~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17a.lmer) ##Age:BMIz p=0.992

##Model 17a - 2/28/17-- ESR1 Site 4 Age Group 2; bmi z-score available for age group 2
AgeGroup17b.lmer<-lmer(ESR1_Site_4~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17b.lmer) ##Age:BMIz p=0.415902

##Model 18a - 2/28/17-- ESR1 Site 1 Age Group 3; bmi z-score available for age group 3
AgeGroup18a.lmer<-lmer(ESR1_Site_1~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18a.lmer) ##Age:BMIz p=0.754

##Model 18b - 2/28/17-- ESR1 Site 4 Age Group 3; bmi z-score available for age group 3
AgeGroup18b.lmer<-lmer(ESR1_Site_4~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18b.lmer) ##Age:BMIz p=0.3153


############### SREBF1 Mean Methylation Models #############
##Model 19 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup19.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup19.lmer) ##Age:BMIz p=0.6870

##Model 20 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup20.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup20.lmer) ##Age:BMIz p=0.73815

##Model 21 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup21.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup21.lmer) ##Age:BMIz p=0.7680

##Model 22 - 2/28/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup22.lmer<-lmer(SREBF1_Site_4~HEI_score+Unmatched_Age*bmiz+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup22.lmer) ##Age:BMIz p=0.7746

##Model 23 - 2/28/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup23.lmer<-lmer(SREBF1_Site_4~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup23.lmer) ##Age:BMIz p=0.783076

##Model 24 - 2/28/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup24.lmer<-lmer(SREBF1_Site_4~HEI_score+Unmatched_Age*bmiz+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup24.lmer) ##Age:BMIz p=0.7680

######### 03/07/17 - Weight Category:Age Interaction term models #########

#Make sure wtstatus_cat has "Normal" as the reference variable
is.factor(BloodMeth1$wtstatus_cat) #TRUE
levels(BloodMeth1$wtstatus_cat)
#[1] "Normal"      "Obese"       "Overweight"  "Underweight"

############### LINE-1 Mean Methylation Models #############
##Model 1 - 03/07/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA available
AgeGroup1.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup1.lmer) ##Age: Obese p-value = 0.2904; Age:overweight p-value = 0.0679; 
#Age:Underweight p-value = 0.1684

##Model 2 - 03/07/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup2.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup2.lmer) ##Age: Obese p-value = 0.004085; Age:overweight p-value = 0.477257; 
#Age:Underweight p-value = 0.868709

##Model 3 - 03/07/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup3.lmer<-lmer(LINE1_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup3.lmer) ##Age: Obese p-value = 0.5773; Age:overweight p-value = 0.5431; 
#Age:Underweight p-value = 0.6136


############### IGF2 Mean Methylation Models #############
##Model 4 - 03/07/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup4.lmer<-lmer(IGF2_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup4.lmer) ##Age: Obese p-value = 0.445547; Age:overweight p-value = 0.552982; 
#Age:Underweight p-value = 0.796548

##Model 5 - 03/07/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup5.lmer<-lmer(IGF2_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup5.lmer) ##Age: Obese p-value = 0.042588; Age:overweight p-value = 0.954196; 
#Age:Underweight p-value = 0.354976

##Model 6 - 03/07/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup6.lmer<-lmer(IGF2_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup6.lmer) ##Age: Obese p-value = 0.2933; Age:overweight p-value = 0.2842; 
#Age:Underweight p-value = 0.8100

############### H19 Mean Methylation Models #############
##Model 7 - 03/07/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup7.lmer<-lmer(H19_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup7.lmer) ##Age: Obese p-value = 0.7810; Age:overweight p-value = 0.7675; 
#Age:Underweight p-value = 0.4121

##Model 8 - 03/07/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup8.lmer<-lmer(H19_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup8.lmer) ##Age: Obese p-value = 0.0757; Age:overweight p-value = 0.01067; 
#Age:Underweight p-value = 0.49816

##Model 9 - 03/07/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup9.lmer<-lmer(H19_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup9.lmer) ##Age: Obese p-value = 0.23960; Age:overweight p-value = 0.49141; 
#Age:Underweight p-value = 0.90069

############### PPARA Mean Methylation Models #############
##Model 10 - 03/07/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup10.lmer<-lmer(PPARA_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup10.lmer) ##Age: Obese p-value = 0.6919; Age:overweight p-value = 0.9728; 
#Age:Underweight p-value = 0.2293

##Model 11 - 03/07/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup11.lmer<-lmer(PPARA_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup11.lmer) ##Age: Obese p-value = 0.1878; Age:overweight p-value = 0.9688; 
#Age:Underweight p-value = 0.0159

##Model 12 - 03/07/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup12.lmer<-lmer(PPARA_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup12.lmer) ##Age: Obese p-value = 0.6670; Age:overweight p-value = 0.5304; 
#Age:Underweight p-value = 0.6652

############### LEP Mean Methylation Models #############
##Model 13 - 03/07/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup13.lmer<-lmer(LEP_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup13.lmer) ##Age: Obese p-value = 0.994138; Age:overweight p-value = 0.743548; 
#Age:Underweight p-value = 0.358250

##Model 14 - 03/07/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup14.lmer<-lmer(LEP_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup14.lmer) ##Age: Obese p-value = 0.6071; Age:overweight p-value = 0.2439; 
#Age:Underweight p-value = 0.3773

##Model 15 - 03/07/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup15.lmer<-lmer(LEP_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup15.lmer) ##Age: Obese p-value = 0.511; Age:overweight p-value = 0.676; 
#Age:Underweight p-value = 0.361


############### ESR1 Mean Methylation Models #############
##Model 16 - 03/07/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##Age: Obese p-value = 0.288; Age:overweight p-value = 0.9727; 
#Age:Underweight p-value = 0.5666

##Model 17 - 03/07/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup17.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17.lmer) ##Age: Obese p-value = 0.4901; Age:overweight p-value = 0.9014; 


##Model 18 - 03/07/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup18.lmer<-lmer(ESR1_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18.lmer) ##Age: Obese p-value = 0.3142; Age:overweight p-value = 0.1137; 
#Age:Underweight p-value = 0.6170

##Break ESR1 down by CpG site due to differences in longitudinal patterns by site.

##Model 16a - 03/07/17 -- ESR1 Site 1 Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Site_1~HEI_score+Unmatched_Age*wtstatus_cat+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##Age: Obese p-value = 0.58378; Age:overweight p-value = 0.48045; 
#Age:Underweight p-value = 0.88843

##Model 16a - 03/07/17 -- ESR1 Site 4 Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup16.lmer<-lmer(ESR1_Site_4~HEI_score+Unmatched_Age*wtstatus_cat+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup16.lmer) ##Age: Obese p-value = 0.7039; Age:overweight p-value = 0.7515; 
#Age:Underweight p-value = 0.3264

##Model 17a - 03/07/17-- ESR1 Site 1 Age Group 2; bmi z-score available for age group 2
AgeGroup17a.lmer<-lmer(ESR1_Site_1~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17a.lmer) ##Age: Obese p-value = 0.8403; Age:overweight p-value = 0.508218; 
#Age:Underweight p-value = N/A

##Model 17a - 03/07/17-- ESR1 Site 4 Age Group 2; bmi z-score available for age group 2
AgeGroup17b.lmer<-lmer(ESR1_Site_4~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup17b.lmer) ##Age: Obese p-value = 0.98410; Age:overweight p-value = 0.38278; 
#Age:Underweight p-value = N/A

##Model 18a - 03/07/17-- ESR1 Site 1 Age Group 3; bmi z-score available for age group 3
AgeGroup18a.lmer<-lmer(ESR1_Site_1~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18a.lmer) ##Age: Obese p-value = 0.167564; Age:overweight p-value = 0.078292; 
#Age:Underweight p-value = 0.463368

##Model 18b - 03/07/17-- ESR1 Site 4 Age Group 3; bmi z-score available for age group 3
AgeGroup18b.lmer<-lmer(ESR1_Site_4~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup18b.lmer) ##Age: Obese p-value = 0.4444; Age:overweight p-value = 0.2592; 
#Age:Underweight p-value = 0.5847


############### SREBF1 Mean Methylation Models #############
##Model 19 - 03/07/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup19.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup19.lmer) ##Age: Obese p-value = 0.4591; Age:overweight p-value = 0.7272; 
#Age:Underweight p-value = N/A

##Model 20 - 03/07/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup20.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup20.lmer) ##Age: Obese p-value = 0.90705; Age:overweight p-value = 0.11961; 
#Age:Underweight p-value = 0.50490

##Model 21 - 03/07/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup21.lmer<-lmer(SREBF1_Mean_Methylation~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup21.lmer) ##Age: Obese p-value = 0.00298; Age:overweight p-value = 0.13813; 
#Age:Underweight p-value = 0.74633

##Model 22 - 03/07/17 -- Age Group 1; wfl z-score available for age group 1; NO MVPA
AgeGroup22.lmer<-lmer(SREBF1_Site_4~HEI_score+Unmatched_Age*wtstatus_cat+childgender+(1|id), data=AgeGroup1) 
summary(AgeGroup22.lmer) ##Age: Obese p-value = 0.8748; Age:overweight p-value = 0.0706; 
#Age:Underweight p-value = N/A

##Model 23 - 03/07/17-- Age Group 2; bmi z-score available for age group 2
AgeGroup23.lmer<-lmer(SREBF1_Site_4~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup2) 
summary(AgeGroup23.lmer) ##Age: Obese p-value = 0.549933; Age:overweight p-value = 0.574475; 
#Age:Underweight p-value = 0.988076

##Model 24 - 03/07/17-- Age Group 3; bmi z-score available for age group 3
AgeGroup24.lmer<-lmer(SREBF1_Site_4~HEI_score+Unmatched_Age*wtstatus_cat+MVPA_perday+childgender+(1|id), data=AgeGroup3) 
summary(AgeGroup24.lmer) ##Age: Obese p-value = 0.00788; Age:overweight p-value = 0.31379; 
#Age:Underweight p-value = 0.94482
