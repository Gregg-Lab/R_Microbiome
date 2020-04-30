#TITLE
###Sac Sheet Plotter

#DESCRIPTION:
###Plots data from Sac Sheets from an Excel file (.xlsx or .xls)
###Format of data sheets MUST be variables as columns (i.e. weight, height, age, etc.) and observations as rows (i.e. mouse #, sample #, human #, etc.)

#AUTHOR:
###Date Created:  4-15-2019
###Date Last Modified: 06-12-2019
###Created by Zach Carlson, BS
###email: zcarlson@umich.edu
###DO NOT REMOVE ABOVE TEXT#
library(readxl) #excel manipulation package
library(lubridate) #convenient date converter
library(tidyverse) #tidys data
library(reprex) #reproducible examples

######READ IN FUNCTIONS FILE######
source("~/Downloads/R_Microbiome/utilities.R")
##################################

setwd("~/Downloads/R_Microbiome/Sac/") #set working directory
data <- read_excel("met_pn.xlsx", skip = 1) #read in excel sheet
data2 <- data[-16, -c(17,18,22)] #row 16 missing, deleting, column 17, 18, 22 missing, deleting
data2$`Ile Mass (mg)` <- as.numeric(data2$`Ile Mass (mg)`)#column is chr, converting to num
plot <- ggplot(data = data2,
               aes(x = Dose, y = `Ins Dose (uL)`, fill = Dose)) +
  geom_boxplot() +
  geom_point()
plot