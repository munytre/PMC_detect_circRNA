library(tidyverse)
library(RVenn)
library(UpSetR)
library(rtracklayer)
library(SummarizedExperiment)
library(ggpubr)

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
                         ~dplyr::select(.x,
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

###Load in custom FA - mapped fastqs -------------------------------------------
#Load in idxstats files
setwd("LOCATION OF CustomFA (MAPPED) IDXSTATS MATE 1")
mate1 <- loadIdx("LOCATION OF CustomFA (MAPPED) IDXSTATS MATE 1")
setwd("LOCATION OF CustomFA (MAPPED) IDXSTATS MATE 2")
mate2 <- loadIdx("LOCATION OF CustomFA (MAPPED) IDXSTATS MATE 2")
#Filter files to get relevant data
CustomFAMp <- filterCustomCountsSE(mate1, mate2) %>%
  map(~subset(.x,
              circ_id %in% circ_names))
rm(mate1, mate2)

###Load in CIRI2 mapped fastqs -------------------------------------------------
#Load in gtf files
setwd("LOCATION OF CIRI2 (PLASMA) GTFs")
CIRI2Mp <- loadGTF("LOCATION OF CIRI2 (PLASMA) GTFs")
CIRI2Mp <- map(CIRI2Mp,
               ~select(.x, circ_id, bsj)) %>%
  map(~subset(.x,
              circ_id %in% circ_names))

###Load in CIRIquant - custom bed ----------------------------------------------
#Load in gtf files
setwd("LOCATION OF CIRIQUANT (PLASMA) GTFs")
CIRIquantMp <- loadGTF("LOCATION OF CIRIQUANT (PLASMA) GTFs")
CIRIquantMp <- map(CIRIquantMp,
                   ~select(.x, circ_id, bsj)) %>%
  map(~subset(.x,
              circ_id %in% circ_names))

###Number comparison -----------------------------------------------------------
#Prepare counts for all 4 methods
quant <- as_tibble(jobs)
colnames(quant) <- "sample"
quant$CIRI2Mp <- map(CIRI2Mp,
                     ~nrow(.x))
quant$CIRIquantMp <- map(CIRIquantMp,
                         ~nrow(.x))
quant$CustomFAMp <- map(CustomFAMp,
                        ~nrow(.x))
quant$CustomFAUnmp <- map(CustomFAUnmp,
                          ~nrow(.x))
quant <- pivot_longer(quant, cols = 2:5, names_to = "method", values_to = "circs")
quant$circs <- as.numeric(quant$circs)
#Plot counts
ggplot(quant, aes(x = sample, y = circs, fill=factor(method))) +
  geom_bar(stat="identity",position="dodge") +
  scale_y_log10() +
  scale_fill_brewer(name = "Method - Bar", palette = "Oranges") +
  stat_summary(fun.y=mean, geom="line", aes(group = method, color = method)) +
  theme_minimal() +
  ylab("# circRNAs") +
  ggtitle("Detection of circRNAs between methods") +
  theme(axis.text.x = element_text(angle=60,
                                   hjust=1)) +
  scale_colour_discrete(name = "Method - Line")
#Plot counts (No barplot)
ggplot(quant, aes(x = sample, y = circs, fill=factor(method))) +
  stat_summary(fun.y=mean, geom="line", aes(group = method, color = method)) +
  theme_minimal() +
  scale_y_log10() +
  ylab("# circRNAs") +
  ggtitle("Detection of circRNAs between methods") +
  theme(axis.text.x = element_text(angle=60,
                                   hjust=1)) +
  scale_colour_discrete(name = "Method")

###UPSET -----------------------------------------------------------------------
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
#Function to isolate all circs in list of data frames
circ_idIsolate <- function(dataset){
  temp <- map(dataset,
              ~pull(.x, circ_id))
  temp <- as_tibble(unique(unlist(temp)))
  return(temp)
}
#Make one big list of all methods with circ_ids
complete_list <- list(CIRI2Mp, CIRIquantMp, CustomFAMp, CustomFAUnmp)
names(complete_list) <- c("CIRI2Mp", "CIRIquantMp", "CustomFAMp", "CustomFAUnmp")
circ_list <- map(complete_list,
                 ~circ_idIsolate(.x))

#Isolate circ_ids per sample
circids <- map(circ_list,
               'value')
#Concatenate all present circ_ids
circs <- concCirc(circids)
#Generate binary table of circ_ids per sample
UpSet_circ <- commBinCircs(circids,
                           circs)
#Plot UpSet
upset(UpSet_circ,
      order.by = c("freq",
                   "degree"),
      decreasing = c(F,
                     T))
#Plot UpSet
upset(UpSet_circ,
      sets = c("CIRI2Mp","CustomFAUnmp"),
      order.by = c("freq",
                   "degree"),
      decreasing = c(F,
                     T))

###COUNT_DISTRIBUTION ----------------------------------------------------------
#Search for overlapping samples in 2 methods
circids_venn <- Venn(circids)
overlap_CIRI2Mp_CustomFAUnmp <- RVenn::overlap(circids_venn, slice = c("CIRI2Mp","CustomFAUnmp"))
overlap_CIRI2Mp <- map(CIRI2Mp,
                       ~filter(.x, circ_id %in% overlap_CIRI2Mp_CustomFAUnmp))
names(overlap_CIRI2Mp) <- jobs
overlap_CustomFAUnmp <- map(CustomFAUnmp,
                            ~filter(.x, circ_id %in% overlap_CIRI2Mp_CustomFAUnmp))
#Filter circRNAs in only diagnostics samples
diag_overlap_CIRI2Mp <- bind_rows(overlap_CIRI2Mp, .id = "sample") %>%
  filter(sample %in% diagnosis_samples)
diag_overlap_CustomFAUnmp <- bind_rows(overlap_CustomFAUnmp, .id = "sample") %>%
  filter(sample %in% diagnosis_samples)
names(diag_overlap_CustomFAUnmp) <- names(diag_overlap_CIRI2Mp)
#Isolate circRNAs with minimum of 2 bsj
diag_overlap_CIRI2Mp <- filter(diag_overlap_CIRI2Mp, bsj >= 2)
diag_overlap_CustomFAUnmp <- filter(diag_overlap_CustomFAUnmp, bsj >= 2)
#Method 1 = CIRI2Mp | Method 2 = CustomFAUnmp
merged_overlap <- bind_rows(diag_overlap_CIRI2Mp, diag_overlap_CustomFAUnmp, .id = "method")
#Plot violin plot
ggplot(merged_overlap, aes(x = method, y = bsj, fill = method)) +
  geom_violin(color = "white") +
  scale_y_log10() +
  theme_minimal() +
  geom_boxplot(width=0.1) +
  theme(legend.position = "none") +
  scale_x_discrete(labels= c("CIRI2Mp","CustomFAUnmp")) +
  labs(y = "log10(#BSJ)", title = "#BSJ in overlapping circRNAs") +
  stat_compare_means(method = "t.test")

###END -------------------------------------------------------------------------