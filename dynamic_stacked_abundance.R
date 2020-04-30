#TITLE
###Dynamic Stacked Community Relative Abundance Plots

#DESCRIPTION:
###Plots community relative abundance plots for multiple taxonomic levels
###Can view either Phylum, Class, Order, and Family.  
###Genus is theoretically possible, but the large legend makes the graph impossible to see.

#AUTHOR:
###Date Created: 06-19-2018
###Date Last Modified:  07-11-2019
###Created by Zach Carlson, BS
###email: zcarlson@umich.edu
###DO NOT REMOVE ABOVE TEXT#

library(RColorBrewer) #gives color schemes
library(tidyverse)
library(gplots) #data plotter
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(tm) #loaded to use removeNumbers() which removes any number in a string
library(shiny)

######ENTER IN VALUES######
TITLE <- "Relative Test"
TAXONOMY_LEVEL <- "Phylum"
TAXONOMY.FILE <- "combined.final.0.03.cons.taxonomy"
TIMEPOINT.FILE <- "combined.final.shared"
DESIGN.FILE <- "combined_P19.design.txt" #Cols: Sample, Group
###########################


######READ IN FUNCTIONS FILE######
source("~/Downloads/R_Microbiome/utilities.R")
##################################

#get into  correct directory
setwd("~/Downloads/R_Microbiome/Stacked_Abundance/") #mac users

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
#add back in OTU column using row names
tax_clean$OTU <- rownames(tax_clean)

###PART 2. CREATE RELATIVE ABUNDANCE FILE###
#selects for samples with 1000+ total reads and otus with 200+ total reads
otu.info <- subset(otu.info, select = -c(label, numOtus))
otu.info <- otu.info[rowSums(otu.info) > 1000, ] 

#select just p19 values
otu.info$Temp_Name <- rownames(otu.info)
otu.info <- otu.info[otu.info$Temp_Name %in% meta$Sample,] #selects the exact number of samples, accounts for otu.info files with multiple timepoints
otu.info <- subset(otu.info, select = -Temp_Name)

#selects correct otus
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
sums.total <- bind_rows(all_SUMS_CALC_TABLES, .id = "group")
sums.total <- as.data.frame(sums.total) #change from matrix

#add back unique sample IDs
group <- sums.total$group 
sums.total <- subset(sums.total, select = -c(group))
ids <- 1:nrow(sums.total)
names <- paste(group, ids)
rownames(sums.total) <- names
sums.total <- t(sums.total)
data.export <- as.data.frame(sums.total)
data.export <- data.export * 100
data.export$Taxa <- rownames(data.export)
write_csv(data.export, path = "sums_total.csv")


#create empty data frame
size <- as.numeric(dim(sums.total)[1] * dim(sums.total)[2])
final.data <- data.frame(matrix(NA, nrow = size, ncol = 3))
colnames(final.data)[1] <- "taxa"
colnames(final.data)[2] <- "names"
colnames(final.data)[3] <- "abundance"

#cycle through relative abundance data and add to new ggplot friendly table
#get dims
GGPLOT.ROWS <- as.numeric(dim(sums.total)[1])
GGPLOT.COLS <- as.numeric(dim(sums.total)[2])

counter = 1 #counter for new data table

for(i in 1:GGPLOT.COLS){
  for(j in 1:GGPLOT.ROWS){
    final.data[counter, "taxa"] <- row.names(sums.total)[j]
    final.data[counter, "names"] <- colnames(sums.total)[i]
    final.data[counter, "abundance"] <- sums.total[j,i]
    counter = counter + 1
  }
}

plots <- ggplot(data=final.data, 
                aes(x=names, y=abundance, fill=taxa)) +
  geom_bar(stat="identity", color="black") +
  labs(title = "Stacked Relative Abundance Plot", y = "Relative Abundances", x = "Groups") +
  coord_flip() +
  theme(legend.position="bottom")
plots
