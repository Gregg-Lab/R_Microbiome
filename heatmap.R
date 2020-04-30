#TITLE
###Heatmap

#DESCRIPTION:
###Plots a section of OTUs in heatmap format

#AUTHOR:
###Date Created:  06-12-2019
###Date Last Modified:  07-09-2019
###Created by Zach Carlson, BS
###email: zcarlson@umich.edu
###DO NOT REMOVE ABOVE TEXT#

#load required libraries
library(RColorBrewer) #gives color schemes
library(tidyverse) #data plotter
library(gplots) #data plotter
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(tm) #loaded to use removeNumbers() which removes any number in a string

######READ IN FUNCTIONS FILE######
source("~/Downloads/R_Microbiome/utilities.R")
##################################

# VARIABLES ##############################################################
TAXONOMY.FILE <- "combined.final.0.03.cons.taxonomy" #Cols: Size, Taxonomy
TAXONOMY_LEVEL <- "Genus"
TIMEPOINT.FILE <- "combined.final.shared" #Cols: Otu0001, Otu0002, etc.
DESIGN.FILE <- "combined_P19.design.txt" #Cols: Sample, Group
TITLE <- "P19 Heatmap"
##########################################################################

#set working directory to where your .shared file is and load file
setwd("~/Downloads/R_Microbiome/Heatmap/")

##import data
tax <- read.table(file=TAXONOMY.FILE, row.names = 1,
                  header=TRUE,
                  check.names=FALSE,
                  comment.char="")
otu.info.raw <- read.table(file=TIMEPOINT.FILE, row.names = 2, header = TRUE)
meta <- read.table(file=DESIGN.FILE, header = TRUE)

#CLEAN TAXONOMY FILE
tax_clean <- cleanTaxonomy(tax, minOtu = 200)
TAX.LIST <- sort(unique(tax_clean[,eval(TAXONOMY_LEVEL)]))
TAX.LENGTH <- as.numeric(length(TAX.LIST))
tax_clean <- subset(tax_clean, select = eval(TAXONOMY_LEVEL))
tax_clean$OTU <- rownames(tax_clean)

###CREATE RELATIVE ABUNDANCE FILE###
otu.info <- subset(otu.info.raw, select = -c(label, numOtus))
otu.info <- otu.info[rowSums(otu.info) > 1000, ] 
otu.info$Temp_Name <- rownames(otu.info)
otu.info <- otu.info[otu.info$Temp_Name %in% meta$Sample,] #selects the exact number of samples, accounts for otu.info files with multiple timepoints
otu.info <- subset(otu.info, select = -Temp_Name)
otu.info <- otu.info[,names(otu.info) %in% tax_clean$OTU]
otu.matrix <- as.matrix(otu.info)
otu.rel <- otu.matrix/rowSums(otu.matrix) #check with 'rowSums(otu.rel)' all rows = 1
otu.rel.max <- apply(otu.rel, 2, max) #2 means columns, apply function max to columns in otu.rel
otu.rel.filtered <- otu.rel[, otu.rel.max>0.02]
otu.rel.filtered <- as.data.frame(otu.rel.filtered)

#add in groups, assign to column names
otu.rel.filtered$Sample <- rownames(otu.rel.filtered)
final.data <- merge(otu.rel.filtered, meta, by = "Sample")
final.data <- final.data[order(final.data$Group),]
names <- as.character(final.data$Group)
final.data <- subset(final.data, select = -c(Sample, Group))
final.data <- as.data.frame(t(final.data))
names(final.data) <- names

#add in tax info, assign to row names
final.data$OTU <- rownames(final.data)
final.data <- merge(final.data, tax_clean, by ="OTU")
genus <- final.data$Genus
final.data <- subset(final.data, select = -c(OTU, Genus))
final.data <- as.matrix(final.data)
final.data <- t(final.data)
colnames(final.data) <- genus
final.data <- final.data[,1:20] #select first fifteen otus
#myCol3 <- colorRampPalette(brewer.pal(8,"GnBu"))(8) #GnBu.. classic green to blue pallette
#myBreaks2 <- c(0, 0.001, 0.003, 0.01, 0.05, 0.10, 0.50, 0.80, 1) 

#heatmap for p19
heatmap.2(final.data, 
          dendrogram = 'none', #removes dendogroms
          key      = TRUE, #removes color key
          #col      = myCol3, #regular heatmap
          #breaks   = myBreaks2,  #regular heatmap
          Colv = FALSE, #doesn't reorder column
          Rowv = FALSE, #doesn't reorder rows so mets and ctrls are separate
          scale    = "column", #used if want red/blue z-distribution
          tracecol = "#303030", #used if want red/blue z-distribution
          col = bluered, #used if want red/blue z-distribution
          main     = "P19 Heatmap", #title
          trace    = 'none',
          cexRow = 1, #sets y axis label text size
          cexCol = 1, #reduces x axis label text size
          margins = c(8,8), #sets graph margins
          srtCol = 30) #adds 30? angle to x axis labels






##LDA Results Specific Heatmap
###############################################################
otu.info.lda <- subset(otu.info.raw, select = c(Otu0007, Otu0002, Otu0006, Otu0010, Otu0011, Otu0005, Otu0023, Otu0021, Otu0022, Otu0017, Otu0018, Otu0031, Otu0027, Otu0026, Otu0033, Otu0041, Otu0024)) 
###############################################################

TAX.LIST <- sort(unique(tax_clean[,eval(TAXONOMY_LEVEL)]))
TAX.LENGTH <- as.numeric(length(TAX.LIST))
rownames(tax_clean) <- tax_clean$OTU
tax_clean <- subset(tax_clean, select = eval(TAXONOMY_LEVEL))
#add back in OTU column using row names
tax_clean$OTU <- rownames(tax_clean)

#CREATE RELATIVE ABUNDANCE FILE###
otu.info.lda <- otu.info.lda[rowSums(otu.info) > 1000, ] 
otu.info.lda$Temp_Name <- rownames(otu.info.lda)
otu.info.lda <- otu.info.lda[otu.info.lda$Temp_Name %in% meta$Sample,] #selects the exact number of samples, accounts for otu.info.lda files with multiple timepoints
otu.info.lda <- subset(otu.info.lda, select = -Temp_Name)
otu.matrix <- as.matrix(otu.info.lda)
otu.rel <- otu.matrix/rowSums(otu.matrix) #check with 'rowSums(otu.rel)' all rows = 1
otu.rel <- as.data.frame(otu.rel)
otu.rel$Sample <- rownames(otu.rel) #add OTUs column so we can merge
final.data <- merge(otu.rel, meta, by = "Sample")
final.data <- final.data[order(final.data$Group),]
names <- as.character(final.data$Group)
final.data <- subset(final.data, select = -c(Sample, Group))
final.data <- as.data.frame(t(final.data))
names(final.data) <- names
final.data$OTU <- rownames(final.data)
final.data <- merge(final.data, tax_clean, by ="OTU", sort = FALSE)
genus <- final.data$Genus
final.data <- subset(final.data, select = -c(OTU, Genus))
final.data <- as.matrix(final.data)
final.data <- t(final.data)
colnames(final.data) <- genus
#myCol3 <- colorRampPalette(brewer.pal(8,"GnBu"))(8) #GnBu.. classic green to blue pallette
#myBreaks2 <- c(0, 0.001, 0.003, 0.01, 0.05, 0.10, 0.50, 0.80, 1) 

#heatmap for p19
heatmap.2(final.data, 
          dendrogram = 'none', #removes dendogroms
          key      = TRUE, #removes color key
          #col      = myCol3, #regular heatmap
          #breaks   = myBreaks2,  #regular heatmap
          Colv = FALSE, #doesn't reorder column
          Rowv = FALSE, #doesn't reorder rows so mets and ctrls are separate
          scale    = "column", #used if want red/blue z-distribution
          tracecol = "#303030", #used if want red/blue z-distribution
          col = bluered, #used if want red/blue z-distribution
          main     = "P19 Heatmap", #title
          trace    = 'none',
          cexRow = 1, #sets y axis label text size
          cexCol = 1, #reduces x axis label text size
          margins = c(8,8), #sets graph margins
          srtCol = 30) #adds 30? angle to x axis labels