library(foreign)
library(stats) 
library(splines) 
library(mgcv) 
library(nlme) 

#setwd("C:/Users/lab/Desktop/Joe/P01 Mouse Stats") ##Set working directory to location of file of interest

setwd("C:/Users/jjkoch/Dropbox/P01 Mouse Tail Stats")

TailMeth1<- read.csv("All_Gene_Tail_Methylation_Avg_022916.csv")

##Version of spreadsheet without mediterranean groups; for DOHaD conference poster.
TailMeth3<- read.csv("All_Gene_Tail_Methylation_Avg_4 Exp Groups_101216.csv")

attach(TailMeth1) ##Set created dataframe as working data set for models
#detach(TailMeth1) ##Remove working dataframe
attach(TailMeth3)

##Create a complete dataset
TailMeth2 <- na.omit(TailMeth3[,c("Exposure","Exposure.Coded","Age", "Gene", "Mean.Methylation")]) 
detach(TailMeth3)
attach(TailMeth2) 
detach(TailMeth2)

library(reshape2)
library(ggplot2)
library(gridExtra)

?ggplot 

##Change order of grouping variable to swap order of facet grid. Otherwise 10M is first.
TailMeth2$Age <- factor(TailMeth2$Age,
                        levels = c("PND21", "10M"))
##Change order of Exposure variable to get Control first in order.
TailMeth2$Exposure <- factor(TailMeth2$Exposure,
                        levels = c("Control", "Control + BPA", "Western", "Western + BPA"))

Esr1 <- TailMeth2[which(TailMeth2$Gene=="Esr1"),]
LINE1 <- TailMeth2[which(TailMeth2$Gene=="LINE1"),]
IAP <- TailMeth2[which(TailMeth2$Gene=="IAP"),]
H19 <- TailMeth2[which(TailMeth2$Gene=="H19"),]
Igf2 <- TailMeth2[which(TailMeth2$Gene=="Igf2"),]

##Create titles with italics
Esr1_title <- expression(paste(italic("Esr1")))
Igf2_title <- expression(paste(italic("Igf2")))
H19_title <- expression(paste(italic("H19")))

##Add to plot to change title (??)
#+ labs(y=Esr1_title)

##Line Plot for Esr1
Esr1_LinePlot<-ggplot(data=Esr1, aes(x=Age, y=Mean.Methylation, group=Exposure, shape=Exposure, color=Exposure)) +
  geom_line(size=0.75) +
  geom_point(aes(shape=Exposure), size=2.5) +  xlab("Age Group") + ylab("Mean Methylation (%)") +
  ggtitle(Esr1_title)+ theme_bw() + theme(legend.position="none") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(colour = "black", size=1.5)) + scale_colour_hue(l=45)+ ylim(0, 16)+
  theme(plot.title = element_text(size=12, face="bold", vjust=2),
        axis.title.x = element_text(size=12, vjust=-0.25),
        axis.title.y = element_text(size=15, vjust=1.25))

##Line Plot for LINE1
LINE1_LinePlot<-ggplot(data=LINE1, aes(x=Age, y=Mean.Methylation, group=Exposure, shape=Exposure, color=Exposure)) +
  geom_line(size=0.75) +
  geom_point(aes(shape=Exposure), size=2.5) +  xlab("Age Group") + ylab("Mean Methylation (%)") +
  ggtitle("LINE-1")+ theme_bw() + theme(legend.position="none") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(colour = "black", size=1.5)) + scale_colour_hue(l=45)+ ylim(63, 65.5)+
  theme(plot.title = element_text(size=12, vjust=2),
        axis.title.x = element_text(size=12, vjust=-0.25),
        axis.title.y = element_text(size=15, vjust=1.25))

##Line Plot for IAP
IAP_LinePlot<-ggplot(data=IAP, aes(x=Age, y=Mean.Methylation, group=Exposure, shape=Exposure, color=Exposure)) +
  geom_line(size=0.75) +
  geom_point(aes(shape=Exposure), size=2.5) +  xlab("Age Group") + ylab("Mean Methylation (%)") +
  ggtitle("IAP")+ theme_bw() + theme(legend.position="none") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(colour = "black", size=1.5)) + scale_colour_hue(l=45)+ ylim(87.5, 91.5)+
  theme(plot.title = element_text(size=12, vjust=2),
        axis.title.x = element_text(size=12, vjust=-0.25),
        axis.title.y = element_text(size=15, vjust=1.25))

##Line Plot for H19
H19_LinePlot<-ggplot(data=H19, aes(x=Age, y=Mean.Methylation, group=Exposure, shape=Exposure, color=Exposure)) +
  geom_line(size=0.75) +
  geom_point(aes(shape=Exposure), size=2.5) +  xlab("Age Group") + ylab("Mean Methylation (%)") +
  ggtitle(H19_title)+ theme_bw() + theme(legend.position="none") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(colour = "black", size=1.5)) + scale_colour_hue(l=45)+ ylim(45, 63)+
  theme(plot.title = element_text(size=12, face="bold", vjust=2),
        axis.title.x = element_text(size=12, vjust=-0.25),
        axis.title.y = element_text(size=15, vjust=1.25))


##Line Plot for Igf2
Igf2_LinePlot<-ggplot(data=Igf2, aes(x=Age, y=Mean.Methylation, group=Exposure, shape=Exposure, color=Exposure)) +
  geom_line(size=0.75) +
  geom_point(aes(shape=Exposure), size=2.5) +  xlab("Age Group") + ylab("Mean Methylation (%)") +
  ggtitle(Igf2_title)+ theme_bw() + theme(legend.position="none") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(colour = "black", size=1.5)) + scale_colour_hue(l=45)+ ylim(30, 40)+
  theme(plot.title = element_text(size=12, face="bold", vjust=2),
        axis.title.x = element_text(size=12, vjust=-0.25),
        axis.title.y = element_text(size=15, vjust=1.25))
Igf2_LinePlot

##Plot all graphs on a single page
grid.arrange(LINE1_LinePlot, IAP_LinePlot, Esr1_LinePlot, H19_LinePlot,Igf2_LinePlot, ncol=3)

##Create separate version of just one graph to extract legend.
Igf2_LinePlot_legend<-ggplot(data=Igf2, aes(x=Age, y=Mean.Methylation, group=Exposure, shape=Exposure, color=Exposure)) +
  geom_line(size=0.75) +
  geom_point(aes(shape=Exposure), size=2.5) +  xlab("Age Group") + ylab("Mean Methylation (%)") +
  ggtitle(Igf2_title)+ theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(colour = "black", size=1.5)) + scale_colour_hue(l=45)+ ylim(30, 40)+
  theme(plot.title = element_text(size=12, face="bold", vjust=2),
        axis.title.x = element_text(size=12, vjust=-0.25),
        axis.title.y = element_text(size=15, vjust=1.25))
Igf2_LinePlot_legend

##Create a common legend for ggplot
library(gridExtra)
##Define function to create legend from plot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

##Step 1: Create plots -- completed above

##Step 2: Define legend
legend <- get_legend(Igf2_LinePlot_legend)
##Increase font size of legend
legend2 <- legend + theme(legend.text=element_text(size=3))

##Step 3: remove legend from all 5 plots
LINE1_LinePlot <- LINE1_LinePlot + theme(legend.position="none")
IAP_LinePlot <- IAP_LinePlot + theme(legend.position="none")
Esr1_LinePlot <- Esr1_LinePlot + theme(legend.position="none")
Igf2_LinePlot <- Igf2_LinePlot + theme(legend.position="none")
H19_LinePlot <- H19_LinePlot + theme(legend.position="none")

##Step 4: Grid arrange all plots with new legend
grid.arrange(LINE1_LinePlot, IAP_LinePlot, Esr1_LinePlot, Igf2_LinePlot, H19_LinePlot, legend, ncol=3)
