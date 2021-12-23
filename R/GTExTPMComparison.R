library(tidyverse)
library(tximportData)
library(tximport)
library(ggpubr)

###RESOURCES -------------------------------------------------------------------
Salmon_outfiles <-  list.files(path = "LOCATION OF SALMON OUTPUT",
                               pattern = ".sf", 
                               full.names = T, 
                               recursive = T)
sample_names <- scan("LOCATION OF GTEX SAMPLE NAMES (RANDOM5)",
                     what = "",
                     quiet=TRUE)
names(Salmon_outfiles) <- sample_names

###FUNCTIONS -------------------------------------------------------------------
#Function to make correlation plots
correlationPlots <- function(output_location, max_samples, select1, select2, samplenames, data){
  #Define output location
  pdf(output_location)
  #Print correlation plots
  for(i in 1:max_samples){
    corplot <- ggplot(data,
                      aes(x=unlist(data[,select1[i]]),
                          y=unlist(data[,select2[i]]))) +
      geom_point() +
      geom_smooth(method = lm,
                  color = "orange") +
      stat_cor(method = "pearson") +
      theme_minimal() +
      labs(title = sample_names[i],
           x = "Salmon",
           y = "GTEx")
    
    print(corplot)
  }
  #Stop printing with dev.off()
  dev.off()
}

###LOAD IN DATA ----------------------------------------------------------------
#Import Salmon quant files, only transcript-level values
txi.salmon.tx <- tximport(Salmon_outfiles,
                          type = "salmon",
                          txOut = T)

#Import the GTEx_TPMs
GTEx_transcript_TPM_Random5 <- read.csv("LOCATION OF GTEX MEDIAN EXPRESSION (RANDOM5)",
                                        sep="")

###PREPARE DATA ----------------------------------------------------------------
#Subset the TPMs from salmon quant and assign new column names
salmon_random5_matrix <- txi.salmon.tx$abundance
salmon_random5 <- as_tibble(txi.salmon.tx$abundance)
colnames(salmon_random5) <- paste0(colnames(salmon_random5),"_salmon")
salmon_random5$transcript_id <- rownames(salmon_random5_matrix)

#Concatenate all TPMs into one data frame
combined_GTEx_salmon_TPMs <- left_join(salmon_random5,
                                       GTEx_transcript_TPM_Random5,
                                       by = "transcript_id") %>%
  drop_na()

###PLOT ------------------------------------------------------------------------
#Make the correlation plots
correlationPlots("OUTPUT LOCATION",
                 5,
                 c(1:5),
                 c(8:12),
                 sample_names,
                 combined_GTEx_salmon_TPMs)

###END -------------------------------------------------------------------------