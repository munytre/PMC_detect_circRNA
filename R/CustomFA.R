library(tidyverse)
library(RVenn)
library(UpSetR)

###RESOURCES -------------------------------------------------------------------
#Jobs file is a file containing sample names of the sequencing data
jobs <- scan("LOCATION OF JOBS FILE (CUSTOMFA)",
             what = "",
             quiet=TRUE)
jobs_with_mates <- scan("LOCATION OF JOBS FILE WITH BOTH SEQUENCING MATES (CUSTOMFA)",
                        what = "",
                        quiet=TRUE)
diagnosis_samples <- c("SAMPLENAME",
                       "SAMPLENAME",
                       "SAMPLENAME",
                       "SAMPLENAME",
                       "SAMPLENAME")

###FUNCTIONS -------------------------------------------------------------------
#Function to load in idxstats-files
loadIdx <- function(path_to_txt){
  #Define list of files
  names <- list.files(path = path_to_txt,
                      pattern = "*.txt")
  
  ##Load in GTFs as one list of tibbles
  txt_list <- map(names,
                  ~read.delim(.,
                              header = F,
                              col.names = c("circ_id","pseudolength", "mapped_reads", "unmapped_reads"),
                              stringsAsFactors = F))
  txt_list <- map(txt_list,
                  as_tibble)
  #Append correct names to the data
  names(txt_list)<-names
  #Return the dataframe
  return(txt_list)
}
#Function to correctly filter sequencing reads from STAR-idxstats (PE-sequencing)
filterCustomCountsPE <- function(Idx_data){
  #Load in idx_data and isolate full read-pairs by division of 2
  custom_counts <- map(Idx_data,
                       ~mutate(.x,
                               mutate(., mapped_reads_unique = mapped_reads/2)))
  #To account for non-paired reads, all values containing 0.5 are floored
  custom_counts <- map(custom_counts,
                       ~mutate(.x,
                               mutate(., reads = floor(mapped_reads_unique))
                       ))
  #Make new list of tibbles containing circ_id and floored read counts
  filtered_counts <- map(custom_counts,
                         ~select(.x,
                                 circ_id, reads))
  #Remove zeroes from the tibbles
  filtered_counts %<>% map(~mutate_at(.x,
                                      vars(reads),
                                      as.integer)) %>%
    map(~na_if(.x, 0)) %>%
    map(~na.omit(.x))
  names(filtered_counts) <- jobs
  #Return list of filtered tibbles
  return(filtered_counts)
}
#Function to correctly filter sequencing reads from STAR-idxstats (SE-sequencing)
filterCustomCountsSE <- function(mate1, mate2){
  joined_counts <- map2(mate1, mate2, inner_join, by = "circ_id")
  full_counts <- map(joined_counts,
                     ~mutate(.x, total_counts = mapped_reads.x + mapped_reads.y))
  #Make new list of tibbles containing circ_id and read counts
  filtered_counts <- map(full_counts,
                         ~select(.x,
                                 circ_id, total_counts))
  #Remove zeroes from the tibbles
  filtered_counts_nonzero <- filtered_counts %>% map(~mutate_at(.x,
                                                        vars(total_counts),
                                                        as.integer)) %>%
    map(~na_if(.x, 0)) %>%
    map(~na.omit(.x))
  
  names(filtered_counts_nonzero) <- jobs
  #Return list of filtered tibbles
  return(filtered_counts_nonzero)
}
#Function to load in concatenated circ_ids
concCirc <- function(circ_ids){
  #Load circ_ids in Venn
  circids_venn <- Venn(circ_ids)
  #Merge unique circ_ids and make from the ids a tibble
  full_list <- RVenn::unite(circids_venn)
  fullset_tibble <- as_tibble(full_list)
  #Return the circ_ids in a tibble
  return(fullset_tibble)
}
#Function to make binary table of common circ_ids
commBinCircs <- function(UpSet_input, concatenated_circ){
  #Load in binary table from each sample
  UpSet_circ <- fromList(UpSet_input)
  #Add full circ_id-list to increase human readability
  UpSet_circ <- as.data.frame(bind_cols(concatenated_circ,
                                        UpSet_circ))
  #Change name of first column to circ_id
  names(UpSet_circ)[1] <- "circ_id"
  #Return complete binary table
  return(UpSet_circ)
}
#Function to plot amount of circRNAs per sample
#reorder = "yes"|||no-reorder = "no"
circFreqPlots <- function(data, plot_title, reorder){
  #Plotting parameters
  if(reorder == "yes"){
    plot <- data %>%
      gather("sample",
             "Frequency",
             -circ_id) %>%
      group_by(sample) %>%
      summarise(n = sum(Frequency)) %>%
      mutate(diagnosis = sample %in% diagnosis_samples) %>%
      ggplot(aes(x = reorder(sample,
                             -n),
                 y = n,
                 fill = diagnosis)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      ggtitle(plot_title) +
      xlab("sample") +
      ylab("# of circRNAs") +
      scale_fill_manual("Diagnosis", values = c("TRUE" = "brown", "FALSE" = "orange")) +
      coord_flip() +
      scale_x_discrete(limits = rev)
    #Return plot
    return(plot)
  }
  else{
    plot <- data %>%
      gather("sample",
             "Frequency",
             -circ_id) %>%
      group_by(sample) %>%
      summarise(n = sum(Frequency)) %>%
      mutate(diagnosis = sample %in% diagnosis_samples) %>%
      ggplot(aes(x = sample,
                 y = n,
                 fill = diagnosis)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      ggtitle(plot_title) +
      xlab("sample") +
      ylab("# of circRNAs") +
      scale_fill_manual("Diagnosis", values = c("TRUE" = "brown", "FALSE" = "orange")) +
      coord_flip() +
      scale_x_discrete(limits = rev)
    #Return plot
    return(plot)
  }
}

###DATA ------------------------------------------------------------------------
#Load in idxstats files
setwd("LOCATION OF CustomFA (UNMAPPED) IDXSTATS MATE 1")
mate1 <- loadIdx("LOCATION OF CustomFA (UNMAPPED) IDXSTATS MATE 1")
setwd("LOCATION OF CustomFA (UNMAPPED) IDXSTATS MATE 2")
mate2 <- loadIdx("LOCATION OF CustomFA (UNMAPPED) IDXSTATS MATE 2")
#Filter files to get relevant data
filtered_counts <- filterCustomCountsSE(mate1, mate2)
#Isolate circ_ids per sample
circids <- map(filtered_counts,
               'circ_id')
#Concatenate all present circ_ids
circs <- concCirc(circids)
#Generate binary table of circ_ids per sample
UpSet_circ <- commBinCircs(circids,
                           circs)

###PLOT ------------------------------------------------------------------------
#Barplot
circFreqPlots(UpSet_circ,
              "#circRNAs per sample (Plasma)",
              "no")
#UpSet-plot
upset(UpSet_circ,
      sets = diagnosis_samples,
      order.by = c("freq",
                   "degree"),
      decreasing = c(F,
                     T))

###SUBSET ----------------------------------------------------------------------
#Subsetting the circ_ids present in all diagnosis samples
circids_venn <- Venn(circids)
diagnosis_shared <- as_tibble(RVenn::overlap(circids_venn,
                                             slice = diagnosis_samples))

###END -------------------------------------------------------------------------