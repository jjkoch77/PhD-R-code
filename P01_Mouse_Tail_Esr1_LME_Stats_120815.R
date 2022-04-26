library(foreign) 
library(stats) 
library(splines) 
library(mgcv) 
library(nlme) 

setwd("C:/Users/lab/Desktop/Joe/P01 Mouse Stats") ##Set working directory to location of file of interest

TailMeth1 <- read.csv("mEsr1_Tail_Methylation_GEE Formatted_120715.csv")

attach(TailMeth1) ##Set created dataframe as working data set for models
#detach(TailMeth1) ##Remove working dataframe


##Create a complete dataset
TailMeth2 <- na.omit(TailMeth1[,c("Litter","id","Sex","Exposure.Group","Exposure_Coded","Age.Group","Age","Mean.Methylation","Site_1", 
                                  "Site_2","Site_3")]) 
detach(TailMeth1)
attach(TailMeth2) 
detach(TailMeth2)

## characteristics at baseline 
summary(Mean.Methylation[Age==1])
summary(Mean.Methylation[Age==0]) ##distribution of methylation at baseline

##Histograms to examine normality
par(mfrow=c(2,1)) ##Makes 2x1 frame histogram. 
hist(Mean.Methylation,nclass=20) ##nclass=20 indicates that you want 20 bars. 

par(mfrow=c(2,2))
hist(Site_1,nclass=20)
hist(Site_2,nclass=20)
hist(Site_3,nclass=20)

par(mfrow=c(1,1))
hist(log(Mean.Methylation),nclass=20) #Checking whether log-transformation helps with normality. In this case, it looks fairly normal before transformation.
##Note: Repeated measures, so log transformation is probably not appropriate here.

##Scatterplots to look at spread of data
plot(Exposure_Coded,Mean.Methylation) ##Scatterplot exposure vs. PND22 methylation 
plot(Exposure_Coded,Site_1) ##Scatterplot exposure vs. Site 1 methylation

boxplot(Mean.Methylation ~ Exposure.Group, main = "Absolute Esr1 Methylation by Exposure Group", 
        ylab = "Absolute Methylation (%)", las = 1)

boxplot(Mean.Methylation ~ Exposure.Group, main = "Absolute Esr1 Methylation by Exposure Group", 
        ylab = "Absolute Methylation (%)", las = 2, cex.axis=0.7, cex=0.75)

##ANOVA analysis between exposure groups at each time point
##First, let's looks at the means across the different groups:
tapply(Mean.Methylation, Exposure_Coded, mean, na.rm=TRUE) #shows means of methylation for different groups

#ANOVa for Mean Methylation across different exposure groups
Exp1<-as.factor(Exposure_Coded)
ANOVA.meth<-aov(Mean.Methylation~Exposure.Group)##ANOVA
summary(ANOVA.meth) 
print(model.tables(ANOVA.meth,"means"),digits=4)       #report the means and the number of subjects/cell
boxplot(Mean.Methylation~Exposure_Coded,data=TailMeth1)
TukeyHSD(ANOVA.meth)


#ANOVa for Site-Specific methylation across different exposure groups
ANOVA.meth1<-aov(Site_1~Exp1)##ANOVA
summary(ANOVA.meth1) 
print(model.tables(ANOVA.meth1,"means"),digits=4)       #report the means and the number of subjects/cell
boxplot(Site_1~Exposure_Coded,data=TailMeth1) 
TukeyHSD(ANOVA.meth1)

ANOVA.meth2<-aov(Site_2~Exposure_Coded)##ANOVA
summary(ANOVA.meth2) 
print(model.tables(ANOVA.meth2,"means"),digits=4)       #report the means and the number of subjects/cell
boxplot(Site_2~Exposure_Coded,data=TailMeth1)
TukeyHSD(ANOVA.meth2)

ANOVA.meth3<-aov(Site_3~Exposure_Coded)##ANOVA
summary(ANOVA.meth3) 
print(model.tables(ANOVA.meth3,"means"),digits=4)       #report the means and the number of subjects/cell
boxplot(Site_3~Exposure_Coded,data=TailMeth1)
TukeyHSD(ANOVA.meth3)

ANOVA.meth4<-aov(Site_4~(Exposure_Coded))##ANOVA
summary(ANOVA.meth4) 
print(model.tables(ANOVA.meth4,"means"),digits=4)       #report the means and the number of subjects/cell
boxplot(Site_4~Exposure_Coded,data=TailMeth1)
TukeyHSD(ANOVA.meth4)

##Create sex factor variable
sex1 <- factor(Sex, labels=c("M", "F"))

######################LME Model - Crude##############################

## Let's do the analysis using lmer from lme4 
install.packages("lme4")
install.packages("lmerTest")
library(lme4) 
library(lmerTest) ##Loading this package adds p-values to summary of lme4 output!! DO THIS.

##Define Control as reference group for Exposure.Group variable. R chose BPA by default.
is.factor(TailMeth2$Exposure_Coded)
levels(TailMeth2$Exposure.Group)
TailMeth2$Exposure.Group1 = factor(TailMeth2$Exposure.Group, c("Control","BPA","Mediterranean", "Mediterranean + BPA", "Western", "Western + BPA"))
levels(TailMeth2$Exposure.Group1)

##Model 1 - ESR1 (12/08/15)
TailMeth.lmer<-lmer(Mean.Methylation~Age+Age*factor(Exposure.Group1)+factor(Sex)+Age*factor(Sex)+(1|id)+(1|Litter), data=TailMeth2) 
summary(TailMeth.lmer)

##Model 2 - ESR1(12/08/15)
TailMeth.lmer1<-lmer(Mean.Methylation~Age+(1|id)+(1|Litter), data=TailMeth2) 
summary(TailMeth.lmer1)

################# --- Exposure Group Separation Models --- ######################

##Creates set of only control samples
names(TailMeth2)
TailMethControl <- (TailMeth2[TailMeth2$Exposure_Coded == 1,])
##Model for control samples ONLY
TailControl.lmer<-lmer(Mean.Methylation~Age+Age*factor(Sex)+(1|id)+(1|Litter), data=TailMethControl) 
summary(TailControl.lmer)

##Scatterplots of data:
detach(TailMeth2)
attach(TailMethControl)
detach(TailMethControl)
plot(Sex, Mean.Methylation) ##Clear difference in methylation levels based on sex.

##Creates set of only Mediterranean samples
names(TailMeth2)
TailMethMed <- (TailMeth2[TailMeth2$Exposure_Coded == 2,])
##Model for control samples ONLY
TailMed.lmer<-lmer(Mean.Methylation~Age+Age*factor(Sex)+(1|id)+(1|Litter), data=TailMethMed) 
summary(TailMed.lmer)

##Creates set of only Med+BPA samples
names(TailMeth2)
TailMethMedBPA <- (TailMeth2[TailMeth2$Exposure_Coded == 3,])
##Model for control samples ONLY
TailMedBPA.lmer<-lmer(Mean.Methylation~Age+Age*factor(Sex)+(1|id)+(1|Litter), data=TailMethMedBPA) 
summary(TailMedBPA.lmer)

##Creates set of only BPA samples
names(TailMeth2)
TailMethBPA <- (TailMeth2[TailMeth2$Exposure_Coded == 4,])
##Model for control samples ONLY
TailBPA.lmer<-lmer(Mean.Methylation~Age+Age*factor(Sex)+(1|id)+(1|Litter), data=TailMethBPA) 
summary(TailBPA.lmer)

##Creates set of only Western samples 
names(TailMeth2)
TailMethWestern <- (TailMeth2[TailMeth2$Exposure_Coded == 5,])
##Model for control samples ONLY
TailWestern.lmer<-lmer(Mean.Methylation~Age+Age*factor(Sex)+(1|id)+(1|Litter), data=TailMethWestern) 
summary(TailWestern.lmer)

##Scatterplots of data:
detach(TailMeth2)
attach(TailMethWestern)
detach(TailMethWestern)

library(reshape2)
library(ggplot2)

##Create boxplot 
Age1<-factor(Age, labels=c("PND22", "10M"))
TailMeth2Plot.Age<-ggplot(TailMethWestern, aes(x = Age1, y = Mean.Methylation)) + geom_boxplot(aes(fill = Sex), width = 0.9) + ggtitle("Absolute Esr1 Methylation on Western Diet") + theme_bw()
TailMeth2Plot.Age

##Add axis labels
TailMeth2Plot.Age + labs(x="Exposure Group", y="Mean All Site Methylation")

##Alter axis labels
TailMeth2Plot.Age + theme(axis.text.x=element_text(angle=90, size=10, vjust=0.5))

##Keep new axis labels, Move axis titles and main title away from plot and change size
TailMeth2Plot.Age + labs(x="Exposure Group", y="Mean All Site Methylation") + theme(plot.title = element_text(size=20, face="bold", vjust=2),
                                                                                    axis.title.x = element_text(vjust=-0.25),
                                                                                    axis.title.y = element_text(vjust=1.25), axis.text.x=element_text(angle=90, size=10, vjust=0.5))

  
##Creates set of only Western+BPA samples
names(TailMeth2)
TailMethWestBPA <- (TailMeth2[TailMeth2$Exposure_Coded == 6,])
##Model for control samples ONLY
TailWestBPA.lmer<-lmer(Mean.Methylation~Age+Age*factor(Sex)+(1|id)+(1|Litter), data=TailMethWestBPA) 
summary(TailWestBPA.lmer)

################ Comparing Models #####################

##ANOVA on model to compare models. 
anova(TailMeth.lmer1, TailMeth.lmer)
##Based on lower AIC, TailMeth.lmer1 -- smaller model -- is a better fit of the data; however, it does not show interaction effects.


##Calculating least squared means and CI for fixed effect factor terms in the model
difflsmeans(TailMeth.lmer1, test.effs=NULL) ##Not sure if useful - 12/08/15



###################Using nlme for mixed effects model. Not as flexible as lme4 in terms of adding random effects.
##Model 1:
Tailmeth.lme<-lme(fixed=Mean.Methylation~Exposure_Coded+Age+sex1, random=~1|id, random=~1|Litter, data=TailMeth2) 
##Model 2:
Tailmeth.lme2<-lme(fixed=Mean.Methylation~Exposure_Coded+Age+sex1, random=~1|Litter, data=TailMeth2) 
##Model 3:
Tailmeth.lme3<-lme(fixed=Mean.Methylation~Exposure_Coded+Age+sex1, random=~1+Litter|id, data=TailMeth2)

summary(Tailmeth.lme) 
summary(Tailmeth.lme2) 
summary(Tailmeth.lme3) 
## diagnostic plots 
plot(sbp.lme)    #default is fitted vs standardized residuals; test for constant variance. 
plot(sbp.lme, locode~resid(.)) ##Residuals by location coding; REsidual = 0.0 at population mean value 
plot(sbp.lme, log(sbp)~fitted(.)|locode, abline=c(0,1)) ##Bottom left = locode 1, top-right = locode 98. 

##Return variance-covariance matrix for random variables 
getVarCov(Tailmeth.lme) ##(Intercept) = variance 

##Extract the coefficients for fixed and random components  
Tailmeth.lme$coef ##Can be used to calculate intercept;EX: For locode 1, intercept = 4.500539 + 5.272532e-03. 
Tailmeth.lme$coef$fixed 
Tailmeth.lme$coef$random 



######################GEE model - crude######################

sex1 <- factor(Sex, labels=c("M", "F"))

PND22_Tail_model1 <-gee(Mean.Methylation~factor(Exposure_Coded)+sex1, id=id, corstr = "independence", data=TailMeth1)

#id=Litter -- cluster by litter number; removes within-litter correlation; useful when looking at single age methylation by exposure group
#id=id -- cluster by sample number; takes into account correlation between repeated measures of methylation.

summary(PND22_Tail_model1)

##Overall GEE model with Age as predictor of methylation.
PND22_Tail_model2 <-gee(Mean.Methylation~Age+sex1+factor(Exposure_Coded), id=id, corstr = "independence", data=TailMeth1)
summary(PND22_Tail_model2)

##Site-specific age predictor GEE models!
PND22_Tail_model3 <-gee(Site_1~Age+sex1+factor(Exposure_Coded), id=id, corstr = "independence", data=TailMeth1)
summary(PND22_Tail_model3)

PND22_Tail_model4 <-gee(Site_2~Age+sex1+factor(Exposure_Coded), id=id, corstr = "independence", data=TailMeth1)
summary(PND22_Tail_model4)

PND22_Tail_model5 <-gee(Site_3~Age+sex1+factor(Exposure_Coded), id=id, corstr = "independence", data=TailMeth1)
summary(PND22_Tail_model5)

PND22_Tail_model6 <-gee(Site_4~Age+sex1+factor(Exposure_Coded), id=id, corstr = "independence", data=TailMeth1)
summary(PND22_Tail_model6)





##GEE Formatting notes from manual
#gee(formula, id,
#data, subset, na.action,
#R = NULL, b = NULL,
#tol = 0.001, maxiter = 25,
#family = gaussian, corstr = "independence",
#Mv = 1, silent = TRUE, contrasts = NULL,
#scale.fix = FALSE, scale.value = 1, v4.4compat = FALSE)