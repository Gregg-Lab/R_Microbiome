#TITLE
###Alpha Diversity Plot

#DESCRIPTION:
###Plots Alpha Diversity metrics

#AUTHOR:
###Date Created:  06-11-2019
###Date Last Modified:  07-11-2019
###Created by Zach Carlson, BS
###email: zcarlson@umich.edu
###DO NOT REMOVE ABOVE TEXT#

#load required libraries
library(RColorBrewer) #gives color schemes
library(tidyverse)
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(tm) #loaded to use removeNumbers() which removes any number in a string
library(shiny)
library(ggpubr)

######READ IN FUNCTIONS FILE######
source("~/Downloads/R_Microbiome/utilities.R")
##################################

######ENTER IN VALUES######
TITLE <- "P11 Mom Alpha Diversity"
SUMMARY.FILE <- "combined_mom.final.p11.groups.summary.csv" #Cols: [standard raw Mothur file]
DESIGN.FILE <- "combined_mom_p11.design.txt" #Cols: Sample, Group
###########################

setwd("~/Downloads/R_Microbiome/AlphaDiversity/")

#import data
stats <- read.table(file = SUMMARY.FILE, header = TRUE)
groups <- read.table(file = DESIGN.FILE, header = TRUE) 

#subset data based on design file samples
groups_ordered<- groups[order(groups$Group, groups$Sample),] #order table for splittin
stats_filtered <- stats[which(stats$group %in% groups$Sample),] #get only samples from stats file that are in design file
groups_filtered <- groups_ordered[which(groups_ordered$Sample %in% stats$group),] #repeat process in reverse, remove all samples from design file that are not in filtered sample (didn't get read, bad sample, etc)

#combine files, select alpha metric, and export
full_data <- merge(groups_filtered, stats_filtered, by.x = "Sample", by.y = "group", sort = FALSE)
shannon_data <- subset(full_data, select = c(Sample, Group, shannon))
write_csv(shannon_data, path = "data_alpha.csv") #export for plotting

#plot alpha metric
ggplot(data = shannon_data,
       mapping = aes(x = Group, y = shannon, fill = Group), mainTitle = TITLE) +
  geom_boxplot() +
  geom_point() 