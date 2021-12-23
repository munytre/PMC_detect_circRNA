library(tidyverse)
library(rtracklayer)
library(SummarizedExperiment)
library(RVenn)
library(UpSetR)
library(gridExtra)

setwd("LOCATION OF CIRI2 (PLASMA) GTFs")

###FUNCTIONS -------------------------------------------------------------------
#Function to load in GTFs
loadGTF <- function(path_to_gtfs){
  #Define list of files
  gtf <- list.files(path = path_to_gtfs,
                    pattern = "*.gtf")
  ##Load in GTFs as one list of tibbles
  gtf_list <- map(gtf,
                  import)
  gtf_list <- map(gtf_list,
                  as_tibble)
  #Append correct names to the data
  names(gtf_list)<-gtf
  #Fix data type of fsj, bsj and junc_ratio
  gtf_list <- map(gtf_list,
                  ~mutate_at(.x,
                             vars(c(fsj,
                                    bsj,
                                    junc_ratio)),
                             as.numeric
                  ))
  #Fix data type of seqnames
  gtf_list <- map(gtf_list,
                  ~mutate_at(.x,
                             vars(seqnames),
                             as.character
                  ))
  #Make list of autosomes and sex chromosomes
  chromosome_selector <- c(1:22)
  chromosome_selector <- append(chromosome_selector,
                                "X")
  #Subset data containing values in chromosome_selector
  gtf_list <- map(gtf_list,
                  ~subset(.x,
                          seqnames %in% chromosome_selector))
  #Return the dataframe
  return(gtf_list)
}
#GENERATION OF ENS_ID ----------------------------------------------------------
#Load in gtfs
gtf_list_full <- loadGTF("LOCATION OF CIRI2 (PLASMA) GTFs")
gtf_list <- gtf_list_full
#Subset circ_ids and ens_id
circids <- map(gtf_list,
               ~select(.x, circ_id, gene_id))
#Concatenate all gtfs
main_list <- circids %>%
  bind_rows() %>%
  unique()
#Filtered circ_ids are loaded in
filtered_circ_ids <- read.table("OUTPUT LOCATION",
                                quote="\"",
                                comment.char="")
#Filter circ_ids
filtered_circ_ENS <- filter(main_list, circ_id %in% filtered_circ_ids$V1)
filtered_ENS <- select(filtered_circ_ENS, gene_id)
#Write table to disk
write.table(filtered_ENS,
      file = "OUTPUT LOCATION",
      col.names = F,
      row.names = F,
      quote = F)

###END -------------------------------------------------------------------------