###TITLE
###Microbiome R Script Template

#DESCRIPTION:
###Description goes here.


#AUTHOR:
###Date Created: MM-DD-YYYY
###Date Last Modified:  MM-DD-YYYY
###[DELETE:  TEMPLATE LAST UPDATED:  06-12-2019]
###Created by Zach Carlson, BS
###email: zcarlson@umich.edu
###DO NOT REMOVE ABOVE TEXT#

###LIBRARIES:
library(RColorBrewer) #gives color schemes
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(tm) #loaded to use removeNumbers() which removes any number in a string
library(shiny)
library(tidyverse) #ggplot2, tibble, purrr, tidyr, readr, dplyr, stringr, forcats
library(ggpubr)

######ENTER IN VALUES######
TITLE <- "Title of Plot"
TAXONOMY_LEVEL <- "Genus" #Kingdom, Phylum, Class, Order, Family, Genus, Species
TAXA_ELEMENT <- "Akkermansia" #Element within a Taxonomy Level
TAXONOMY.FILE <- "file.taxonomy" #Cols: Size, Taxonomy
TIMEPOINT.FILE <- "file.shared" #Cols: OTU0001, OTU0002, etc.
DESIGN.FILE <- "file.txt" #Cols: Row, Sample, Group
###########################

######READ IN FUNCTIONS FILE######
source("~/Downloads/R_Microbiome/utilities.R")
##################################