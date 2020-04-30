#TITLE
###PCoA Plot

#DESCRIPTION:
###Plots a PCoA plot separating by experimental groups

#AUTHOR:
###Date Created: 06-04-2019
###Date Last Modified:  07-09-2019
###Created by Zach Carlson, BS
###email: zcarlson@umich.edu
###DO NOT REMOVE ABOVE TEXT#

library(RColorBrewer) #gives color schemes
library(gplots) #data plotter
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(tm) #loaded to use removeNumbers() which removes any number in a string
library(shiny)
library(tidyverse) #ggplot2, tibble, purrr, tidyr, readr, dplyr, stringr, forcats
library(ggpubr)

######READ IN FUNCTIONS FILE######
source("~/Downloads/R_Microbiome/utilities.R")
##################################

# VARIABLES #######################################################
TAXONOMY.FILE <- "combined_pn.final.0.03.cons.taxonomy" #Cols: Size, Taxonomy
CORR.FILE <- "combined_pn.final.p19.pearson.corr.axes" #Cols: OTU, axis1, p-value, axis2, p-value, length
PCOA.FILE <- "combined_pn.final.0.03.pick.thetayc.0.03.lt.pcoa.axes" #Cols: group, axis1, axis2, etc.
DESIGN.FILE <- "combined_P19.design.txt" #Cols: Sample, Group
TIMEPOINT.FILE <- "combined_pn.final.shared" #Cols: Otu0001, Otu0002, etc.
TITLE <- "P19 PCoA"
###################################################################

setwd("~/Downloads/R_Microbiome/PCoA/")

##import data
tax <- read.table(file=TAXONOMY.FILE, row.names = 1,
                  header=TRUE,
                  check.names=FALSE,
                  comment.char="")
shared <- read.table(file=TIMEPOINT.FILE, row.names = 2, header = TRUE)
design <- read.table(file=DESIGN.FILE, header = TRUE)

#clean data
tax_clean <- cleanTaxonomy(tax)
tax_clean$OTU <- rownames(tax_clean)
shared <- cleanShared(shared, design, tax_clean)
design <- cleanDesign(design, shared)
corr.axes <- read.table(file=CORR.FILE, header = TRUE)
corr.tax <- merge(corr.axes, tax_clean, by = "OTU", all.y = TRUE)
corr.tax.1 <- subset(corr.tax, p.value < 0.001)
dim(corr.tax.1) #25 OTUs
corr.tax.2 <- subset(corr.tax.1, p.value.1 < 0.001)
dim(corr.tax.2) #0 OTUs

##make plot
pcoa.axes <- read.table(file=PCOA.FILE, header = TRUE)
pcoa.pch <- c(1, 16) #1 is open circle, 16 is filled circle
pcoa.pre.final <- merge(design, pcoa.axes, by.x = "Sample", by.y = "group")
pcoa.final <- pcoa.pre.final[order(pcoa.pre.final$Group),]
shape <- data.frame(c(rep(1, 10), rep(16, 8)))  #controls as open, and mets as filled
colnames(shape) <- "shape" #maintains column name of "shape"
pcoa.final.shape <- cbind(shape, pcoa.final) #combines shape labels with data table
#write.csv(pcoa.axes.final, file = "pcoa.axes.final.csv") #exports csv file to run in mothur!

#p19 plot (ALL DATA KEPT)
plot(pcoa.final.shape$axis1, 
     pcoa.final.shape$axis2, 
     pch = pcoa.final.shape$shape, #selects shape from shape column
     col = "black",  #selects color from group column
     main = "P19 PCoA", 
     xlab ="Axis 1 (34.44%)", 
     ylab = "Axis 2 (29.13%)", 
     ylim = c(-.8, .8),  #sets limits from -60% to 60%
     xlim = c(-.8, .8), #sets limits from -60% to 60%
     cex = 1.3)

#add legend
legend(x = -0.3, #places legend and (0.33, .4)
       y = .7, 
       legend = c("Ctrl PN", "Met PN"), 
       pch = pcoa.pch,
       col = c("black"),
       cex = 0.8)
