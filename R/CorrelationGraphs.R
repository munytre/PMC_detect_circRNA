library(tidyverse)
library(ggpubr)
library(DESeq2)

###FUNCTIONS -------------------------------------------------------------------
countToMatrix <- function(location, sample_name){
  #Read in the table
  fcount_table <- read.table(
    location,
    header=T)
  #Assign $Geneid as rownames
  row.names(fcount_table) <- fcount_table$Geneid
  #Remove extra columns and rename them accordingly with sample_name
  fcount_table <- fcount_table[, -c(1:6)]
  fcount_table_names <- as.list(
    t(
      read.table(
        sample_name
      )))
  colnames(fcount_table) <- fcount_table_names
  #Return the dataframe
  return(fcount_table)
}
countDualMData <- function(sample_name, condition1, condition2, n_samples){
  #Load in samplenames
  metadata_table1 <- read.table(
    sample_name,
    col.names = "sample"
  )
  metadata_table2 <- read.table(
    sample_name,
    col.names = "sample"
  )
  #Assign condition to each sample
  metadata_table1$condition <- rep((condition1),
                                   each=n_samples)
  metadata_table2$condition <- rep((condition2),
                                   each=n_samples)
  #Generate full metadata table
  metadata_table <- rbind(metadata_table1,
                          metadata_table2)
  #Return the dataframe
  return(metadata_table)
}
totalCounts <- function(matrix1, matrix2){
  #Merging by rows
  totalCounts <- merge(matrix1,
                       matrix2,
                       by=0,
                       all=T)
  #Convert first column to rownames
  totalCounts <- column_to_rownames(totalCounts,
                                    var="Row.names")
  #Return the dataframe
  return(totalCounts)
}
countNormalization <- function(counts_combined, metadata){
  #Load in data
  dds <- DESeqDataSetFromMatrix(countData = counts_combined,
                                colData = metadata,
                                design = ~ condition)
  #Generate size factors for normalization
  dds <- estimateSizeFactors(dds)
  #Performing median of ratios method of normalization (DESeq2)
  normalized_counts <- counts(dds,
                              normalized=T)
  #Return the dataframe
  return(normalized_counts)
}
MRMCombine <- function(HISAT2, STAR, metadata){
  combined_counts <- totalCounts(HISAT2,STAR)
  MRM_counts <- as.data.frame(countNormalization(combined_counts, metadata))
}
correlationPlots <- function(output_location, max_samples, select1, select2, samplenames, MRM){
  #Define output location
  pdf(output_location)
  #Print correlation plots
  for(i in 1:max_samples){
    selector1 <- select1
    selector2 <- select2
    sample_names <- as.list(t(read.table(
      samplenames
    )))
    corplot <- MRM %>%
      ggplot(aes(x=MRM[,selector1[i]], y=MRM[,selector2[i]])) +
      geom_point() +
      geom_smooth(method = lm, color = "orange") +
      stat_cor(method = "pearson") +
      theme_minimal() +
      labs(title = sample_names[i], x = "HISAT2", y = "STAR")
    print(corplot)
  }
  #Stop printing with dev.off()
  dev.off()
}

###DATA ------------------------------------------------------------------------
#Load in counts for bwa and STAR
HISAT2_counts <- countToMatrix(
  "LOCATION OF HISAT2 COUNTS",
  "LOCATION OF JOBS FILE")
STAR_counts <- countToMatrix(
  "LOCATION OF STAR COUNTS",
  "LOCATION OF JOBS FILE")
#Combine the metadata of bwa and STAR
MData <- countDualMData(
  "LOCATION OF JOBS FILE",
  "HISAT2",
  "STAR",
  38
)

#PLOTS -------------------------------------------------------------------------
#Combine normalized counts and return the object as a data frame
MRM_counts <- MRMCombine(HISAT2_counts,
                         STAR_counts,
                         MData)
#Make correlation plots
correlationPlots("OUTPUT LOCATION",
                 38,
                 c(1:38),
                 c(39:76),
                 "LOCATION OF JOBS FILE",
                 MRM_counts)

###END -------------------------------------------------------------------------