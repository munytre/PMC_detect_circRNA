library(tidyverse)
library(rtracklayer)
library(SummarizedExperiment)
library(RVenn)
library(UpSetR)
library(gridExtra)
library(GenomicRanges)
library(ggbeeswarm)
library(eulerr)

###FUNCTIONS -------------------------------------------------------------------
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
  SharedFrequency <- dplyr::select(SharedFrequency,
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
#Function to plot amount of circRNAs per sample
circFreqPlots <- function(gtf_set, plot_title){
  plot <- gtf_set %>%
    gather("sample",
           "Frequency",
           -circ_id) %>%
    group_by(sample) %>%
    summarise(n = sum(Frequency)) %>%
    ggplot(aes(x = reorder(sample,
                           -n),
               y = n)) +
    geom_col(fill = "orange") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
                                                         ends = "last"))) +
    # theme(axis.text.x = element_text(angle=45,
    # hjust=1)) +
    ggtitle(plot_title) +
    xlab("Patient samples (Plasma)") +
    ylab("# of circRNAs")
  # geom_text(aes(label=n), vjust = -1, size=2.75)
  return(plot)
}

###CHROMOSOME SELECTOR ---------------------------------------------------------
#Make list of autosomes and sex chromosomes
chromosome_selector <- c(1:22)
chromosome_selector <- append(chromosome_selector,
                              "X")

###Load in CIRI2 mapped fastqs (Neuroblastoma) ---------------------------------
#Load in GTFs and make tibble of found circ_ids
setwd("LOCATION OF CIRI2 (PRIMARY) GTFs")
NBs <- loadGTF("LOCATION OF CIRI2 (PRIMARY) GTFs")
#Subset data containing values in chromosome_selector
gtf_list_NB <- map(NBs,
                   ~subset(.x,
                           seqnames %in% chromosome_selector))
gtf_list_NB <- map(gtf_list_NB,
                   ~subset(.x,
                           circ_type %in% "exon"))
circids_NB <- map(gtf_list_NB,
                  'circ_id')
genes_NB <- map(gtf_list_NB,
                'gene_id')

###Load in CIRI2 mapped fastqs (Plasma) ----------------------------------------
#Load in GTFs and make tibble of found circ_ids
setwd("LOCATION OF CIRI2 (PLASMA) GTFs")
Plasma <- loadGTF("LOCATION OF CIRI2 (PLASMA) GTFs")
#Subset data containing values in chromosome_selector
gtf_list_plasma <- map(Plasma,
                       ~subset(.x,
                               seqnames %in% chromosome_selector))
gtf_list_plasma <- map(gtf_list_plasma,
                       ~subset(.x,
                               circ_type %in% "exon"))
circids_plasma <- map(gtf_list_plasma,
                      'circ_id')

###MERGING COUNTS --------------------------------------------------------------
#Merge n of circRNAs per sample type
lengths_NB <- as_tibble(lengths(circids_NB))
colnames(lengths_NB) <- "circRNAs"
lengths_Plasma <- as_tibble(lengths(circids_plasma))
colnames(lengths_Plasma) <- "circRNAs"
joined_lengths <- bind_rows(lengths_NB,
                            lengths_Plasma,
                            .id = "source")
joined_lengths[joined_lengths == 1] <- "Neuroblastoma"
joined_lengths[joined_lengths == 2] <- "Plasma"
NBs_n <- joined_lengths[1:51,]
Plasma_n <- joined_lengths[52:75,]
###PLOT ------------------------------------------------------------------------
#Save plot NB as a grid object
grid1 <- ggplot(NBs_n,
                aes(source,
                    circRNAs)) +
  geom_boxplot(fill = "#BA313D") +
  theme_minimal() +
  ylab("# circRNAs") +
  ylim(0,30000) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
#Save plot Plasma as a grid object
grid2 <- ggplot(Plasma_n,
                aes(source,
                    circRNAs)) +
  geom_boxplot(fill = "orange") +
  theme_minimal() +
  ylab("# circRNAs") +
  ylim(0,30000) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
#Plot both grid objects
grid.arrange(grid1,
             grid2,
             ncol=2,
             top = 'CircRNAs detected by CIRI2')

###UNIQUE ----------------------------------------------------------------------
#Isolate unique circRNAs
NB_genes_unique <- unique(unlist(genes_NB))
NB_genes_unique <- unique(unlist(strsplit(NB_genes_unique,",")))
#Load in circRNAs of interest
plasma_circs <- scan("LOCATION OF PLASMA CIRCRNAS",
                     what = "",
                     quiet=TRUE)
neuroblastoma_circs <- read.delim("LOCATION OF FULL CONSENSUS SUBSET") %>% pull(circ_id)
#Find intersection between unique circRNA sets
NB_unique <- unique(unlist(circids_NB))
Plasma_unique <- unique(unlist(circids_plasma))
shared_NB_Plasma <- intersect(NB_unique,
                    Plasma_unique)
#Find intersection between circRNAs of interest and unique sets
shared_filter <- intersect(neuroblastoma_circs,
                           plasma_circs)
#Set up euler (Venn) object
euler_shared <- euler(c(Neuroblastoma = length(NB_unique) - length(shared_NB_Plasma),
                        Plasma = length(Plasma_unique) - length(shared_NB_Plasma),
                        "Neuroblastoma&Plasma" = length(shared_NB_Plasma)),
                      shape = "ellipse")
#Plot eulers (Venn)
plot(euler_shared,
     fills = c("#aa5866","#ffb366"),
     quantities = c(length(NB_unique) - length(shared_NB_Plasma),
                    length(Plasma_unique) - length(shared_NB_Plasma),
                    length(shared_NB_Plasma)))
#Set up euler (Venn) object (After filtering)
euler_shared_filt <- euler(c(Neuroblastoma = length(neuroblastoma_circs) - 5598,
                        Plasma = 0,
                        "Neuroblastoma&Plasma" = 5598),
                      shape = "ellipse")
#Plot eulers (Venn)
plot(euler_shared_filt,
     fills = c("#aa5866","#ffb366"),
     quantities = c(length(neuroblastoma_circs) - 5598,
                    0,
                    5598))

###END -------------------------------------------------------------------------