library(tidyverse)
library(rtracklayer)
library(SummarizedExperiment)
library(RVenn)
library(UpSetR)
library(ggpubr)

setwd("LOCATION OF CIRI2 (PRIMARY) GTFs")

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
  #Return the dataframe
  return(gtf_list)
}
#Function to load in concatenated circ_ids
concCirc <- function(circ_ids){
  circids_venn <- Venn(circ_ids)
  full_list <- RVenn::unite(circids_venn)
  fullset_tibble <- as_tibble(full_list)
  return(fullset_tibble)
}
#Function to make binary list of common circ_ids
commBinCircs <- function(UpSet_input, concatenated_circ){
  UpSet_circ <- fromList(UpSet_input)
  UpSet_circ <- as.data.frame(bind_cols(concatenated_circ,
                                        UpSet_circ))
  names(UpSet_circ)[1] <- "circ_id"
  return(UpSet_circ)
}
#Function to generate frequency of circRNAs
SharedFreq <- function(binary_table){
  SharedFrequency <- UpSet_circ
  SharedFrequency$freq <- rowSums(SharedFrequency[-1] == 1)
  SharedFrequency <- select(SharedFrequency,
                            circ_id,
                            freq)
  return(SharedFrequency)
}

###RAW DATA --------------------------------------------------------------------
#Load in GTFs and make tibble of found circ_ids
gtf_list_full <- loadGTF("LOCATION OF CIRI2 (PRIMARY) GTFs")

###CHROMOSOME_SELECTOR ---------------------------------------------------------
#Make list of autosomes and sex chromosomes
chromosome_selector <- c(1:22)
chromosome_selector <- append(chromosome_selector,
                              "X")
#Subset data containing values in chromosome_selector
gtf_list <- map(gtf_list_full,
                ~subset(.x,
                        seqnames %in% chromosome_selector))
###EXON_SELECTOR ---------------------------------------------------------------
#Subset exons from data
gtf_list <- map(gtf_list,
                ~subset(.x,
                        circ_type %in% "exon"))
###SHARED_SELECTOR -------------------------------------------------------------
#Subset all circ_ids from data
circids <- map(gtf_list,
               'circ_id')
#Make tibble of concatenated circ_ids
concatenated_circs <- concCirc(circids)
#Make binary list of common circ_ids
UpSet_circ <- commBinCircs(circids,
                           concatenated_circs)
#Generate shared frequency table
SharedFrequency <- SharedFreq(UpSet_circ)

###MINBSJ_SELECTOR -------------------------------------------------------------
gtf_list <- map(gtf_list,
                ~filter(.x,
                        bsj >= 2))

###FULL LIST OF CIRC_IDs
#Copy concatenated_circs for later table-join
full_circ <- concatenated_circs
#Change "value" into "circ_id"
colnames(full_circ) <- "circ_id"

###FULL LIST OF BSJs -----------------------------------------------------------
#Iterate through the list of gtfs and selecting only bsj's and their circ_ids
#to be merged with full list of circ_ids. After merging, the list of data frames
#is reduced to one data frame.
gtf_bsj_joined <- map(gtf_list,
               ~select(.x,
                       circ_id,
                       bsj) %>%
               left_join(full_circ,
                         .x,
                         by = "circ_id")) %>%
  reduce(full_join,
         by = "circ_id")
#After reduce, colnames needs to be re-appended
colnames(gtf_bsj_joined)[2:ncol(gtf_bsj_joined)] <- names(gtf_list)
#Calculation of the means of bsj per circ_id, excluding NA-values
gtf_bsj_joined$mean <- rowMeans(gtf_bsj_joined[-1],
                                na.rm = T)

###FULL LIST OF CPMs -----------------------------------------------------------
#Iterate through the list of gtfs and selecting only bsj's and their circ_ids
#to be merged with full list of circ_ids. After merging, the list of data frames
#is reduced to one data frame.
gtf_score_joined <- map(gtf_list,
                      ~select(.x,
                              circ_id,
                              score) %>%
                        left_join(full_circ,
                                  .x,
                                  by = "circ_id")) %>%
  reduce(full_join,
         by = "circ_id")
#After reduce, colnames needs to be re-appended
colnames(gtf_score_joined)[2:ncol(gtf_score_joined)] <- names(gtf_list)
#Calculation of the means of bsj per circ_id, excluding NA-values
gtf_score_joined$mean <- rowMeans(gtf_score_joined[-1],
                                  na.rm = T)

###LIST OF CIRC_IDs WITH ENSEMBL-IDs OF HOST-GENE & FREQUENCY IN SAMPLES -------
#Selecting only circ_ids and gene_ids for merge.
gtf_list_sub <- map(gtf_list,
                    ~select(.x,
                            circ_id,
                            gene_id))
#Join all circ_ids and gene_ids, removing NAs and duplicated values
gtf_list_geneIDs <- map(gtf_list_sub,
                        ~left_join(full_circ,
                                   .x,
                                   by = "circ_id")) %>%
  bind_rows() %>%
  na.omit() %>%
  unique()
#Creates column indicating if circ_id has multiple host-genes
gtf_list_geneIDs <- mutate(gtf_list_geneIDs,
                           multiple_ENS = case_when(grepl(",",
                                                          gene_id) ~ 1,
                                                    !grepl(",", gene_id) ~ 0))
#Merge frequency per sample per circ_id with main data frame
gtf_list_geneIDs <- left_join(gtf_list_geneIDs,
                              SharedFrequency,
                              by = "circ_id")

###FINAL DATA ------------------------------------------------------------------
#Merging circ_ids, ENSEMBL_IDs, Sample frequency, bsj's, and mean of bsj's.
BSJ_table <- left_join(gtf_list_geneIDs,
                       gtf_bsj_joined,
                       by = "circ_id")
SCORE_table <- left_join(gtf_list_geneIDs,
                         gtf_score_joined,
                         by = "circ_id")

###OUTPUT ----------------------------------------------------------------------
#Write tables to disk
write.table(BSJ_table,
            file = "OUTPUT LOCATION",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
write.table(SCORE_table,
            file = "OUTPUT LOCATION",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

###CORRELATION TABLE -----------------------------------------------------------
#Prepare data frame with only bsj and score
bsj_score <- as_tibble(BSJ_table$mean)
bsj_score$score <- SCORE_table$mean
#Plot
ggplot(bsj_score,
       aes(x = value,
           y = score)) +
  geom_point() +
  geom_smooth(method="lm",
              se=TRUE,
              fullrange=FALSE,
              level=0.95) +
  stat_cor(method = "pearson") +
  labs(x = "bsj",
       y = "score (CPM)") +
  theme_minimal()

###END -------------------------------------------------------------------------