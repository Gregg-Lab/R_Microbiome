###TITLE
###Microbiome Utility Functions

#DESCRIPTION:
###R file dedicated to commonly used functions when processing microbiome data.

#AUTHOR:
###Date Created: 05-28-2019
###Last Modified:  07-02-2019
###Created by Zach Carlson, BS
###email: zcarlson@umich.edu
###DO NOT REMOVE ABOVE TEXT#

###LIBRARIES:
library(RColorBrewer) #gives color schemes
library(gplots) #data plotter
library(vegan) #ecology package, diversity analysis, etc.
library(plyr) #tool for splitting, applying, and combining data
library(tm) #loaded to use removeNumbers() which removes any number in a string
library(shiny)
library(tidyverse) #ggplot2, tibble, purrr, tidyr, readr, dplyr, stringr, forcats
library(ggpubr)

#Returns a subset of data.frame x selecting for the taxa element y
calculateAbundance <- function(df,taxa_element){
  df <- as.data.frame(df)
  df <- df[eval(taxa_element),]
  return(df)
}


#returns table with just ratio values of Bacteroidetes and Firmicutes bacteria
calculateBFRatioTable <- function(list){
  y <- as.data.frame(list)
  original_nrow <- nrow(y)
  original_ncol <- ncol(y)
  y$Bacteroidetes_Sum <- rep(0, times = nrow(list)) #create GP sums table
  y$Firmicutes_Sum <- rep(0, times = nrow(list)) #create GN sums table
  y$Ratio <- rep(0, times = nrow(list))
  
  #sort by GP and GN
  for(i in 1:original_nrow){
    for(j in 1:original_ncol){
      if(colnames(y[j]) == "Bacteroidetes"){
        y[i, "Bacteroidetes_Sum"] <- y[i,j] + y[i, "Bacteroidetes_Sum"]
      } else if(colnames(y[j]) == "Firmicutes") {
        y[i, "Firmicutes_Sum"] <- y[i,j] + y[i, "Firmicutes_Sum"]
        
      } 
    }
  }
  
  #do Ratio calculation
  for(i in 1:original_nrow){
    y[i, "Ratio"] <- y[i, "Bacteroidetes_Sum"] / y[i, "Firmicutes_Sum"]
  }
  
  return(y)
}


#returns table with just ratio values of GP and GN bacteria
calculateGPGNRatioTable <- function(list){
  y <- as.data.frame(list)
  original_nrow <- nrow(y)
  original_ncol <- ncol(y)
  y$GP <- rep(0, times = nrow(list)) #create GP sums table
  y$GN <- rep(0, times = nrow(list)) #create GN sums table
  y$Ratio <- rep(0, times = nrow(list))
  
  #sort by GP and GN
  for(i in 1:original_nrow){
    for(j in 1:original_ncol){
      if(colnames(y[j]) == "Actinobacteria" | colnames(y[j]) == "Firmicutes"){
        y[i, "GP"] <- y[i,j] + y[i, "GP"]
      } else {
        y[i, "GN"] <- y[i,j] + y[i, "GN"]
      }
    }
  }
  
  #do Ratio calculation
  for(i in 1:original_nrow){
    y[i, "Ratio"] <- y[i, "GP"] / y[i, "GN"]
  }
  
  return(y)
}


#returns sums table with calculated values
calculateSumsTable <- function(list){
  original_length <- ncol(otu.matrix) #num of OK OTUs
  list <- as.data.frame(list)
  names(list)[1:original_length] <- rm_g[,eval(TAXONOMY_LEVEL)]
  for(i in 1:nrow(list)){
    for(j in 1:original_length){
      for(k in 1:length(TAX.LIST)){
        if(colnames(list)[j] == TAX.LIST[k]){
          list[i,original_length+k] <- list[i,j] + list[i,original_length+k]
        }
      }
    }
  }
  
  return(list)
}


#REQUIRES:  -User-created design file, x, with columns "Sample", "Group
#           -R-edited .shared file, y, created from cleanShared
#MODIFIES:  meta data table
#EFFECTS:   -Removes samples that have under minSample total sequence
#              reads, default is 1000.
cleanDesign <- function(design, shared){
  shared$Temp_Name <- rownames(shared) #adds temp column to check meta with
  design <- design[design$Sample %in% shared$Temp_Name,] #removes all files from design file that do not have 1000 reads
  return(design)
}


#REQUIRES:  -Unedited, raw .shared file, x.
#           -User-created design file, y, with columns "Sample", "Group
#           -R-edited .taxonomy file, z, created from cleanTaxonomy, 
#              with at least "OTU" column and one taxonomic level 
#              column
#MODIFIES:  otu.info data table
#EFFECTS: -Removes unnecessary columns "label" and "numOtus"
#         -Removes samples that have under minSample total sequence 
#              reads, default is 1000.
#         -Removes samples that are not in the meta data table
cleanShared <- function(shared, design, taxonomy, minSample = 1000){
  
   shared <- subset(shared, select = -c(label, numOtus)) #deletes label, numOtus columns
   ##ORIGINAL POSITION FOR SUBSETTING 1000##
   shared$Temp_Name <- rownames(shared) #adds temp column to check meta with
   design <- design[design$Sample %in% shared$Temp_Name,] #removes all files from design file that do not have 1000 reads
   shared <- shared[shared$Temp_Name %in% design$Sample,] #removes all shared samples that are not in the design file, useful for timepoint-specific inquiries
   shared <- subset(shared, select = -Temp_Name) #deletes Temp_Name column
   shared <- shared[,names(shared) %in% taxonomy$OTU] #deletes all OTUs with < 200 reads
   shared <- shared[rowSums(shared) > 1000,] #remove all samples without 1000 reads
   
   return(shared)
}


#REQUIRES:  .taxonomy file with "Size" and "Taxonomy" columns.
#MODIFIES:  taxonomy data table
#EFFECTS: -Parses raw .taxonomy file into its taxonomic levels
#         -Removes OTUs that have under minOtu total sequence reads, default is 200
#         -Keeps the pattern "_unclassified"! Use cleanTaxonomyKeepUnclassified() in combination with cleanUnclassified()to remove all OTUs with an unclassified read result --Depreciated as of 7.2.19
#NOTE: minSamples > minOtu.  Bad samples do not contribue to minOtu count.
cleanTaxonomy <- function(taxonomy, 
                          minOtu = 200){
  
  taxonomy <- taxonomy[taxonomy$Size > minOtu,] #selects otus with over 200 hits
  taxonomy <- separate(taxonomy, Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
  taxonomy <- subset(taxonomy, select = -Size) #remove size
  
  for(i in 1:nrow(taxonomy)){
    for(j in 1:ncol(taxonomy)){
      taxonomy[i,j] <- removeNumbers(taxonomy[i,j]) #remove numbers
      taxonomy[i,j] <- str_replace_all(taxonomy[i,j], "[[:punct:]]", "") #removes "(" ")" and all special characters
      taxonomy[i,j] <- str_replace_all(taxonomy[i,j], pattern = " ", replacement = "") #remove extra white space
    }
  }
  return(taxonomy)
}


# #REQUIRES:  .taxonomy file with "Size" and "Taxonomy" columns.
# #MODIFIES:  taxonomy data table
# #EFFECTS: -Parses raw .taxonomy file into its taxonomic levels
# #         -Removes OTUs that have under minOtu total sequence reads, default is 200
# #         -DOESNT remove the pattern "unclassified"  ***Use cleanUnclassified() afterwards
# #NOTE: minSamples > minOtu.  Bad samples do not contribue to minOtu count.
# cleanTaxonomyKeepUnclassified <- function(taxonomy, minOtu = 200){
#   #remove unclassified rows
#   taxonomy <- taxonomy[taxonomy$Size > minOtu,] #selects otus with over 200 hits
#   taxonomy <- separate(taxonomy, Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
#   taxonomy <- subset(taxonomy, select = -Size) #remove size
#   
#   for(i in 1:nrow(taxonomy)){
#     for(j in 1:ncol(taxonomy)){
#       taxonomy[i,j] <- removeNumbers(taxonomy[i,j]) #remove numbers
#       taxonomy[i,j] <- str_replace_all(taxonomy[i,j], "[[:punct:]]", " ") #removes "(" ")" and all special characters
#       taxonomy[i,j] <- str_replace_all(taxonomy[i,j], pattern = " ", replacement = "") #remove extra white space
#     }
#   }
#   taxonomy$OTU <- rownames(taxonomy)
#   return(taxonomy)
# }
# 
# 
# #REQUIRES:  cleaned .taxonomy file, from cleanTaxonomyKeepUnclassified()
# #MODIFIES:  taxonomy data table
# #EFFECTS:  -Removes all rows which have "unclassified" at a specified level (i.e. column)
# #         
# #NOTE: minSamples > minOtu.  Bad samples do not contribue to minOtu count.
# cleanUnclassified <- function(taxonomy, column = "Phylum"){
#   taxonomy <- taxonomy %>% filter(!str_detect(taxonomy[,eval(column)], 'unclassified'))
#   
#   
#   return(taxonomy)
# }


#Creates empty sums table based on TAXONOMY_LEVEL selected filled with sums columns, one for each element in TAX.LIST. 
createEmptySumsTable <- function(list){
  y <- subset(list, select = -c(Sample, Group))
  y <- as.matrix(y)
  class(y) <- "numeric" #change matrix to numeric form rather than character
  colnames(y) <- rm_g[,eval(TAXONOMY_LEVEL)] #removes numbers from colnames. e.g. Bacteroidetes.1 -> Bacteroidetes
  y <- as.data.frame(y)
  for(j in 1:length(TAX.LIST)){
    col.name <- TAX.LIST[j]
    y$temp <- rep.int(0, nrow(y))
    names(y)[ncol(y)] <- col.name
  }
  return(y)
}


#returns the sums table exclusively
exportSumsTable <- function(list){
  original_length <- ncol(otu.matrix)
  y <- list[,(original_length+1):(original_length+TAX.LENGTH)] #adds 1 to get the indelist of the FIRST new column
  return(y)
}

