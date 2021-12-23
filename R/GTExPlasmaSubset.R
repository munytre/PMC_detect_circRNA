library(tidyverse)
library(RVenn)
library(UpSetR)
library(rtracklayer)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(Cairo)
library(magick)
library(circlize)
library(viridis)
library(ggbeeswarm)

###RESOURCES -------------------------------------------------------------------
#Jobs file is a file containing sample names of the sequencing data
jobs <- scan("LOCATION OF JOBS FILE",
             what = "",
             quiet=TRUE)
diagnosis_samples <- c("SAMPLENAME",
                       "SAMPLENAME",
                       "SAMPLENAME",
                       "SAMPLENAME",
                       "SAMPLENAME")
circ_names <- read.delim("LOCATION OF CIRCRNA_IDS FROM CONSENSUS SUBSET",
                         header=T) %>%
  pull(circ_id)
circ_ENS <- read.delim("LOCATION OF HOSTGENES_ENS_IDS FROM CONSENSUS SUBSET",
                       header=T)
NB_subset <- read.table("LOCATION OF CIRCRNA FROM ENRICHED_PLASMA SUBSET",
                        header=FALSE,
                        col.names = c("circ_id",
                                      "gene_id"))
numerical_names <- scan("LOCATION OF CIRCRNA FROM ENRICHED_NUMERICAL SUBSET",
                        what = "",
                        quiet=TRUE)
numerical_subset <- read.csv("LOCATION OF CIRCRNA (FULL DATA) FROM ENRICHED_NUMERICAL SUBSET",
                            sep="",
                            stringsAsFactors=FALSE)
DE_genes <- scan("LOCATION OF DE GENES",
                 what = "",
                 quiet=TRUE)
DE_genes_NB_RESTRICTED_SUBSET <- scan("LOCATION OF DE GENES IN PLASMA",
                                      what = "",
                                      quiet=TRUE)
full_set <- read.delim("LOCATION OF FULL ENRICHED_PLASMA SUBSET")

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
#Function to generate frequency of circRNAs
SharedFreq <- function(binary_table){
  SharedFrequency <- UpSet_circ
  SharedFrequency$freq <- rowSums(SharedFrequency[-1] == 1)
  SharedFrequency <- select(SharedFrequency,
                            circ_id,
                            freq)
  return(SharedFrequency)
}
#Function to generate plot of SharedFreq-results
freqTablePlot <- function(freqTable){
  plot <- ggplot(freqTable,
                 aes(freq)) + 
    geom_bar(fill = "orange") +
    theme_minimal() +
    ylab("# circRNAs") +
    xlab("# samples") +
    ggtitle("Common circRNAs") +
    geom_text(stat='count', aes(label=..count..), vjust = -1, size=2.75)
  return(plot)
}
subsetOverlap <- function(circ_ids, selected_sample){
  circids_venn <- Venn(circ_ids)
  samples_overlap <- overlap(circids_venn,slice = selected_sample)
  samples <- CIRI2_low_expr_hclust[selected_sample]
  overlapping <- map(samples,
                     ~filter(.x,
                             circ_id %in% samples_overlap))
  return(samples_overlap)
}

##Load in custom FA - unmapped fastqs -----------------------------------------
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
###Load in CIRI2 mapped fastqs (PRIMARY) ---------------------------------------
#Load in gtf files
setwd("LOCATION OF CIRI2 (PRIMARY) GTFs")
gtf_list_full <- loadGTF("LOCATION OF CIRI2 (PRIMARY) GTFs")
CIRI2_PRIM_filtered <- map(gtf_list_full,
                           ~select(.x, circ_id, bsj, gene_id, score)) %>%
  map(~subset(.x,
              circ_id %in% circ_names))
###ANALYSIS --------------------------------------------------------------------
plasma_circs <- c(CustomFAUnmp_circids,CIRI2_circids) %>%
  unique()

Circs_found_in_plasma <- subset(NB_subset,
                                circ_id %in% plasma_circs)

NB_found_in_plasma <- subset(Circs_found_in_plasma,
                             gene_id %in% numerical_names)

name_NB_found_in_plasma <- as.character(pull(NB_found_in_plasma,
                                             circ_id))

NB_found_in_plasma <- NB_found_in_plasma %>%
  mutate(CIRI2 = circ_id %in% CIRI2_circids) %>%
  mutate(CustomFasta = circ_id %in% CustomFAUnmp_circids)

gene_interest <- DE_genes

gene_interest_circRNAs <- subset(NB_found_in_plasma,
                                 gene_id %in% gene_interest)

table(gene_interest_circRNAs$gene_id)

###BSJ COUNTS ------------------------------------------------------------------
CIRI2Mp_df <- bind_rows(CIRI2_NB_filtered,
                        .id = "sample")
CIRI2Mp_df <- CIRI2Mp_df[,1:3]

CustomFA_df <- bind_rows(CustomFAUnmp,
                         .id = "sample")
colnames(CustomFA_df) <- c("sample",
                           "circ_id",
                           "bsj")

total_bsj_df <- bind_rows(CIRI2Mp_df,
                          CustomFA_df,
                          .id = "method")
total_bsj_df$method <- sub("2", "CustomFA_Plasma", total_bsj_df$method)
total_bsj_df$method <- sub("1", "CIRI2_Plasma", total_bsj_df$method)

sample_annot <- read.delim("LOCATION OF METADATA FOR PLASMA")

total_bsj_df$SampleId <- sub("-.*", "", total_bsj_df$sample)

total_bsj_df <- left_join(total_bsj_df,
                          sample_annot[,c(1,5)],
                          by = "SampleId")
total_bsj_df <- total_bsj_df[,c(1,2,3,4,6)]

total_bsj_df$sample <- sub("-.*", "", total_bsj_df$sample)

CIRI2PRIM_df <- bind_rows(CIRI2_PRIM_filtered,
                          .id = "sample")
CIRI2PRIM_df <- CIRI2PRIM_df[,1:3]
CIRI2PRIM_df$TimepointId <- 0
CIRI2PRIM_df$method <- "CIRI2_Primary"

CIRI2PRIM_df$sample <- sub("\\.gtf", "", CIRI2PRIM_df$sample)

PRIM_PLAS_df <- bind_rows(total_bsj_df,
                          CIRI2PRIM_df)

PRIM_PLAS_DE_circs <- filter(PRIM_PLAS_df,
                       circ_id %in% gene_interest_circRNAs$circ_id)
PRIM_PLAS_DE_circs <- left_join(PRIM_PLAS_DE_circs,
                                full_set[,c(1,2,5)],
                                by = "circ_id")

ggplot(data = PRIM_PLAS_DE_circs,
       aes(TimepointId,
           bsj,
           color = circ_id)) +
  geom_point(aes(shape = method)) +
  scale_y_log10() +
  theme_minimal() +
  facet_wrap(~ gene_name,
             ncol = 3) +
  labs(title = "Detected BSJs in Neuroblastoma and Plasma RNA sequencing data",
       y = "BSJs",
       x = "Timepoint",
       color = "IDs of circRNA",
       shape = "Detection method")
  

### ----------------------------------------------------------------------------
filtered_total_bsj_df <- filter(total_bsj_df,
                                circ_id %in% gene_interest_circRNAs$circ_id)

filtered_total_bsj_df <- left_join(filtered_total_bsj_df,
                                   circ_ENS,
                                   by = "circ_id")

avg_bsj_df <- filtered_total_bsj_df %>%
  group_by(circ_id) %>%
  summarise(avg = mean(bsj))

avg_bsj_df <- left_join(avg_bsj_df,
                        circ_ENS,
                        by = "circ_id")

ggplot(filtered_total_bsj_df,
       aes(circ_id, bsj)) +
  geom_boxplot()

###NB_RESTRICTED FOUND IN PLASMA -----------------------------------------------
NB_RESTRICTED_found_in_plasma <- subset(Circs_found_in_plasma,
                                        gene_id %in% DE_genes_NB_RESTRICTED_SUBSET)

unique(NB_RESTRICTED_found_in_plasma$gene_id)

intersect(unique(gene_interest_circRNAs$gene_id), unique(NB_RESTRICTED_found_in_plasma$gene_id))

write.table(NB_RESTRICTED_found_in_plasma,
            file = "OUTPUT LOCATION",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
###END -------------------------------------------------------------------------