#TITLE
###Dynamic Relative Abundance Plots

#DESCRIPTION:
###Plots relative abundance plots for multiple taxonomic levels
###First, specify either Phylum, Class, Order, and Family.  
###Second, specify which taxa to plot.  
###Finally, Adds stats test
###E.g. LEVEL: Phylum; TAXA: Firmicutes


#AUTHOR:
###Date Created: 4-23-2019
###Last Modified: 07-11-2019
###Created by Zach Carlson, BS
###email: zcarlson@umich.edu
###DO NOT REMOVE ABOVE TEXT#

library(RColorBrewer) #gives color schemes
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(tm) #loaded to use removeNumbers() which removes any number in a string
library(shiny)
library(tidyverse) #ggplot2, tibble, purrr, tidyr, readr, dplyr, stringr, forcats
library(ggpubr)

######READ IN FUNCTIONS FILE######
source("~/Downloads/R_Microbiome/utilities.R")
##################################

######ENTER IN VALUES######
TITLE <- "Dynamic Abundance Plot"
TAXONOMY_LEVEL <- "Class"
TAXA_ELEMENT <- "Bacilli"
TAXONOMY.FILE <- "cohort_d.final.0.03.cons.taxonomy" #Cols: Size, Taxonomy
TIMEPOINT.FILE <- "cohort_d.final.shared" #Cols: OTU0001, OTU0002, etc.
DESIGN.FILE <- "cohort_d_4wk.design.txt" #Cols: Sample, Group
###########################

#get into correct directory
setwd("~/Downloads/R_Microbiome/Individual_Abundance/") #mac users

#READ IN RAW FILES
tax <- read.table(file=TAXONOMY.FILE, 
                  row.names = 1,
                  header=TRUE,
                  check.names=FALSE,
                  comment.char="") #will become otu.class
otu.info <- read.table(file=TIMEPOINT.FILE, header=TRUE, row.names = 2) #will become rm_g
meta <- read.table(file=DESIGN.FILE, row.names = NULL, header =TRUE)

#CLEAN TAXONOMY FILE
tax_clean <- cleanTaxonomy(tax, minOtu = 200)

###GET TAX LEVEL AND THE CORRESPONDING ELEMENTS ###
TAX.LIST <- sort(unique(tax_clean[,eval(TAXONOMY_LEVEL)]))
TAX.LENGTH <- as.numeric(length(TAX.LIST))
###############

#SUBSET TAXONOMY
tax_clean <- subset(tax_clean, select = eval(TAXONOMY_LEVEL))
tax_clean$OTU <- rownames(tax_clean)


###PART 2. CREATE RELATIVE ABUNDANCE FILE###
###PART 2A.  CLEAN FILES (remove bad OTUs and samples)
#selects for samples with 1000+ total reads and otus with 200+ total reads
otu.info <- subset(otu.info, select = -c(label, numOtus))
otu.info <- otu.info[rowSums(otu.info) > 1000, ] #get rid of bad samples
#get a list of "good" samples
otu.info$Temp_Name <- rownames(otu.info)

#removes any sample in meta that didn't have enough reads in otu.info
meta <- meta[meta$Sample %in% otu.info$Temp_Name,] 

#removes all samples in otu.info that are not in the specific meta file (in other words, if the otu.info has multiple timepoints, meta can be used to subset by timepoint)
otu.info <- otu.info[otu.info$Temp_Name %in% meta$Sample,]
otu.info <- subset(otu.info, select = -Temp_Name)

#finally, delete all OTUs that were removed in master tax_clean file
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

#export .csv files of all sum tables
for(i in names(all_SUMS_CALC_TABLES)){
  write_csv(all_SUMS_CALC_TABLES[[i]], path = paste0(i,".csv"))
}

#take all sums tables and combine for plotting and exporting
sums.total <- bind_rows(all_SUMS_CALC_TABLES, .id="group") #returns the name of the list item each observation belonged to, in this case, group.
sums.total <- as.data.frame(sums.total) #change from matrix
data.export <- sums.total
data.export[,2:ncol(data.export)] <- data.export[,2:ncol(data.export)] * 100 #multiplies everything by 100 to convert to percentage.  excludes first column (is group column with character value)
write_csv(data.export, path = "sums_total.csv")



#CALCULATE TAXA ELEMENT AND MAKE GGPLOT COMPATIBLE TABLE
#DOES ONE AT A TIME BASED ON VALUE OF TAXA_ELEMENT AND TAXONOMY_LEVEL
selected_taxa_group <- subset(sums.total, select = group)
selected_taxa_element <- subset(sums.total, select = eval(TAXA_ELEMENT))
selected_taxa_full <- cbind(selected_taxa_group, selected_taxa_element)
selected_taxa_full$taxa <- eval(TAXA_ELEMENT)
names(selected_taxa_full)[2] <- "abundance"
selected_taxa_full$abundance <- selected_taxa_full$abundance * 100 #convert to percentage

plots <- ggplot(data=selected_taxa_full, 
                mapping = aes(x=group, y=abundance, fill=group)) +
  geom_boxplot() +
  #facet_wrap(~taxa) +
  labs(title = "4wk c_Bacilli Relative Abundance Plot", y = "Relative Abundances (%)", x = "Groups") + 
  
  stat_compare_means(comparisons = list(c("Ctrl", "HFD", "Met")))

plots

