#TITLE
###GP:GN Ratio

#DESCRIPTION:
###Group phyla based on presence of Gram-Positive or Gram-Negative bacteria.  
###Assuming both are mutually exclusive.  
###Gram-Positive:  Actinobacteria, Firmicutes
###Gram-Negative:  All other phyla.

#SOURCE:  
###PMID:  17400789

#AUTHOR:
###Date Created:  4-15-2019
###Date Last Modified: 07-11-2019
###Created by Zach Carlson, BS
###email: zcarlson@umich.edu
###DO NOT REMOVE ABOVE TEXT#

library(RColorBrewer) #gives color schemes
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(tm) #loaded to use removeNumbers() which removes any number in a string
library(shiny)
library(tidyverse) #ggplot2, tibble, purrr, tidyr, readr, dplyr, stringr, forcats

######READ IN FUNCTIONS FILE######
source("~/Downloads/R_Microbiome/utilities.R")
##################################

######ENTER IN VALUES######
TITLE <- "GP:GN Ratio"
TAXONOMY_LEVEL <- "Phylum"
TAXONOMY.FILE <- "combined.final.0.03.cons.taxonomy"
TIMEPOINT.FILE <- "combined.final.shared"
DESIGN.FILE <- "combined_8WK.design.txt" #Cols: Sample, Group
###########################

#get into correct directory
setwd("~/Downloads/R_Microbiome/GPGNRatio/") #mac users

#read in raw files 
tax <- read.table(file=TAXONOMY.FILE, 
                  row.names = 1,
                  header=TRUE,
                  check.names=FALSE,
                  comment.char="") #will become otu.class
otu.info <- read.table(file=TIMEPOINT.FILE, header=TRUE, row.names = 2) #will become rm_g
meta <- read.table(file=DESIGN.FILE, row.names = NULL, header =TRUE)


#CLEAN TAXONOMY FILE
tax_clean <- cleanTaxonomy(tax, minOtu = 200) 

###GET PHYLUM###
TAX.LIST <- sort(unique(tax_clean$Phylum))
TAX.LENGTH <- as.numeric(length(TAX.LIST))

#SUBSET TAXONOMY
tax_clean <- subset(tax_clean, select = eval(TAXONOMY_LEVEL))
#add back in OTU column using row names
tax_clean$OTU <- rownames(tax_clean)

###PART 2. CREATE RELATIVE ABUNDANCE FILE###
otu.info <- subset(otu.info, select = -c(label, numOtus))

#select just 8WK values
otu.info$Temp_Name <- rownames(otu.info)
otu.info <- otu.info[otu.info$Temp_Name %in% meta$Sample,]
otu.info <- subset(otu.info, select = -Temp_Name)

#selects samples with over 1000+ total sequence reads
otu.info <- otu.info[rowSums(otu.info) > 1000, ] 



#selects correct otus (> 200 total reads across samples)
otu.info <- otu.info[,names(otu.info) %in% tax_clean$OTU]

otu.matrix <- as.matrix(otu.info)
otu.rel <- otu.matrix/rowSums(otu.matrix) #check with 'rowSums(otu.rel)' all rows = 1
#otu.rel <- subset(otu.rel, select = -c(label, numOtus))
otu.rel.t <- t(otu.rel) #transpose #check with 'colSums(otu.rel.t)' all columns = 1
otu.rel.t <- as.data.frame(otu.rel.t) #make it so we can add back in OTU names without changing to list
otu.rel.t$OTU <- rownames(otu.rel.t) #add OTUs column so we can merge
rm_g <- merge(tax_clean, otu.rel.t, by.x = "OTU", by.y = "OTU")
###PART 2. COMPLETE###
### rm_g FILE COMPLETED ###


otubar <- as.matrix(rm_g)
names(otubar) <- names(rm_g) #as.matrix deletes col.names and row.names, immediately reassign column names because no ordering has occured.
rownames(otubar) <- rm_g[,eval(TAXONOMY_LEVEL)]
otubar <- subset(otubar, select = -c(1:2)) #deletes first two columns.  Col1 = OTU, Col2 = eval(TAXONOMY_LEVEL) selection
bar <- as.data.frame(t(otubar))
bar$Sample <- rownames(bar)
bar_meta <- merge(meta, bar, by.x = "Sample", by.y = "Sample")
bar_ordered <- bar_meta[order(bar_meta$Group),] 

#splits into separate groups
##GET GROUPS##
GROUP.LIST <- as.character(sort(unique(meta$Group)))
GROUP.LENGTH <- as.numeric(length(GROUP.LIST))
###

all <- split(bar_ordered, bar_ordered$Group) 
#list2env(all, envir=.GlobalEnv) #converts lists to data frames

#create sums columns in all tables
all_SUMS <- lapply(all, createEmptySumsTable) #function defined at top of script

#calculate sums columns in all tables
all_SUMS_CALC <- lapply(all_SUMS, calculateSumsTable) #fxn defined at top of script

#export sums tables
all_SUMS_CALC_TABLES <- lapply(all_SUMS_CALC, exportSumsTable) #fxn defined at top of script

all_SUMS_CALC_TABLES_RATIO <- lapply(all_SUMS_CALC_TABLES, calculateGPGNRatioTable) #fxn defined at top of script

#export .csv files of all sum tables
for(i in names(all_SUMS_CALC_TABLES_RATIO)){
  write_csv(all_SUMS_CALC_TABLES_RATIO[[i]], path = paste0(i,".csv"))
}

#take all sums tables and combine for plotting
sums.total <- bind_rows(all_SUMS_CALC_TABLES_RATIO, .id = "group")
sums.total <- as.data.frame(sums.total) #change from matrix
sums.total <- sums.total[,c("GP", "GN", "Ratio", "group")]

#FINAL EXPORT CSV
write.csv(sums.total, "GP_GNRATIO.csv")

ggplot(data = sums.total,
       mapping = aes(x = group, y = Ratio)) +
  geom_boxplot()