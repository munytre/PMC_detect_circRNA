library(tidyverse)
library(RVenn)
library(UpSetR)
library(rtracklayer)
library(SummarizedExperiment)

###RESOURCES -------------------------------------------------------------------
#Jobs file is a file containing sample names of the sequencing data
jobs <- scan("LOCATION OF JOBS FILE (CUSTOMFA)",
             what = "",
             quiet=TRUE)
diagnosis_samples <- c("SAMPLENAME",
                       "SAMPLENAME",
                       "SAMPLENAME",
                       "SAMPLENAME",
                       "SAMPLENAME")
circ_names <- read.delim("LOCATION OF CIRCRNA_IDS FROM CONSENSUS SUBSET",
                         header=FALSE) %>%
  pull(V1)

Gent_count_table <- read.delim("LOCATION OF BIOFLUID ATLAS (PLATELET POOR PLASMA)")

Gent_reformatted_names <- read.delim("LOCATION OF FULL BIOFLUID ATLAS",
                                     header = F)
hg38_circBase_circs <- read.table("LOCATION OF CIRCBASE CIRCRNAS",
                                  quote="\"",
                                  comment.char="")
hg38_circAtlas_circRNAs <- read.table("LOCATION OF CIRCATLAS CIRCRNAS",
                                      quote="\"",
                                      comment.char="")
TSCD_fetal_circRNA <- read.table("LOCATION OF TSCD FETAL CIRCRNAS",
                                 quote="\"",
                                 comment.char="")
TSCD_adult_circRNA <- read.table("LOCATION OF TSCD ADULT CIRCRNAS",
                                 quote="\"",
                                 comment.char="")
hg38_circRNADb_circRNA <- read.table("LOCATION OF CIRCRNADB CIRCRNAS",
                                     quote="\"",
                                     comment.char="")
hg38_circBank_circRNA <- read.table("LOCATION OF CIRCBANK CIRCRNAS",
                                    quote="\"",
                                    comment.char="")
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
  #Load circ_ids in Venn
  circids_venn <- Venn(circ_ids)
  #Merge unique circ_ids and make from the ids a tibble
  full_list <- RVenn::unite(circids_venn)
  fullset_tibble <- as_tibble(full_list)
  #Return the circ_ids in a tibble
  return(fullset_tibble)
}

###Load in custom FA - unmapped fastqs -----------------------------------------
#Load in idxstats files
setwd("LOCATION OF CustomFA (UNMAPPED) IDXSTATS MATE 1")
mate1 <- loadIdx("LOCATION OF CustomFA (UNMAPPED) IDXSTATS MATE 1")
setwd("LOCATION OF CustomFA (UNMAPPED) IDXSTATS MATE 2")
mate2 <- loadIdx("LOCATION OF CustomFA (UNMAPPED) IDXSTATS MATE 2")
#Filter files to get relevant data
CustomFAUnmp <- filterCustomCountsSE(mate1, mate2) %>%
  map(~subset(.x,
              circ_id %in% circ_names))
rm(mate1, mate2)

CustomFAUnmp_circids <- map(CustomFAUnmp,
                            'circ_id') %>%
  concCirc() %>%
  pull()

###Load in CIRI2 mapped fastqs -------------------------------------------------
#Load in gtf files
setwd("LOCATION OF CIRI2 (PLASMA) GTFs")
CIRI2Mp <- loadGTF("LOCATION OF CIRI2 (PLASMA) GTFs")
names(CIRI2Mp) <- jobs
CIRI2_NB_filtered <- map(CIRI2Mp,
                         ~select(.x, circ_id, bsj, gene_id, score)) %>%
  map(~subset(.x,
              circ_id %in% circ_names))

CIRI2_circids <- map(CIRI2_NB_filtered,
                     'circ_id') %>%
  concCirc() %>%
  pull()

###CONCATENATE CIRCS -----------------------------------------------------------
plasma_circs <- c(CustomFAUnmp_circids,CIRI2_circids) %>%
  unique()

###PREP COUNT TABLE ------------------------------------------------------------
Gent_table <- Gent_count_table
Gent_table$nonBED_ID <- Gent_reformatted_names$V1

count_table_NA_ex <- filter(Gent_table,
                            PPP_1 + PPP_2 != 0)

count_table_circs <- pull(count_table_NA_ex,
                          var = 'nonBED_ID')

###CHECK OVERLAP ---------------------------------------------------------------
test <- plasma_circs %in% count_table_circs
summary(test)
intersected_circs <- intersect(plasma_circs,
                               count_table_circs)
intersected_circs

###CIRCS OF INTEREST OVERLAP ---------------------------------------------------
DE_genes_circs <- read.delim("LOAD IN DE GENES WITH CIRCRNA IN PLASMA")
#Subset of genes with found circRNAs
found_circs <- c(unique(DE_genes_circs$circ_id))
#Overlapping circs
failed_circs <- intersect(intersected_circs,
                          found_circs)
#For all circs
failed_circs_full_biofluid <- intersect(found_circs,
                            c(Gent_reformatted_names$V1))

#Comparison with circBase
failed_circs_circBase <- intersect(found_circs,
                                   c(hg38_circBase_circs$V1))

#Comparison with circAtlas
failed_circs_circAtlas <- intersect(found_circs,
                                   c(hg38_circAtlas_circRNAs$V1))
#Comparison with TSCD_fetal
failed_circs_TSCD_fetal <- intersect(found_circs,
                                     c(TSCD_fetal_circRNA$V1))
#Comparison with TSCD_adult
failed_circs_TSCD_adult <- intersect(found_circs,
                                     c(TSCD_adult_circRNA$V1))
#Comparison with circRNADb
failed_circs_circRNADb <- intersect(found_circs,
                                    c(hg38_circRNADb_circRNA$V1))
#Comparison with circBank
failed_circs_circBank <- intersect(found_circs,
                                    c(hg38_circBank_circRNA$V1))
#Generation of presence table
presence_table <- as_tibble(DE_genes_circs[,1:2])
presence_table <- presence_table %>%
  mutate(biofluid_atlas_PPP = circ_id %in% failed_circs) %>%
  mutate(biofluid_atlas = circ_id %in% failed_circs_full_biofluid) %>%
  mutate(circBase = circ_id %in% failed_circs_circBase) %>%
  mutate(circAtlas = circ_id %in% failed_circs_circAtlas) %>%
  mutate(circRNADb = circ_id %in% failed_circs_circRNADb) %>%
  mutate(circBank = circ_id %in% failed_circs_circBank) %>%
  mutate(TSCD_fetal = circ_id %in% failed_circs_TSCD_fetal) %>%
  mutate(TSCD_adult = circ_id %in% failed_circs_TSCD_adult)

###SAVE PLASMA_CIRCS -----------------------------------------------------------
# write.table(plasma_circs,
#             file = "OUTPUT LOCATION",
#             quote = F,
#             row.names = F,
#             col.names = F)

###SAVE Presence_table ---------------------------------------------------------
write.table(presence_table,
            file = "OUTPUT LOCATION",
            quote = F,
            row.names = F,
            col.names = T)

###END -------------------------------------------------------------------------