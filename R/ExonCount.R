library(tidyverse)

###PLASMA CIRCRNA ID -----------------------------------------------------------
#Load in ID of circRNAs found in plasma
circ_plasma <- read.delim("LOCATION OF CIRCRNA_ID",
                          header=FALSE) %>%
  pull()

###EXON COUNTS -----------------------------------------------------------------
#Load in exon count data
exon_count <- read.table("LOCATION OF EXON COUNT DATA", quote="\"", comment.char="")
count <- as_tibble(table(exon_count))
colnames(count) <- c("exon_count", "freq")
count$exon_count <- str_replace(count$exon_count, ">", "")
#Isolate exon numbers for plasma circRNAs
count_plasma <- count %>%
  filter(exon_count %in% circ_plasma)
#Plot for NB-set
ggplot(count, aes(x = freq)) +
  geom_bar(fill = "orange") +
  theme_minimal() +
  ylab("freq") +
  xlab("# exons") +
  ggtitle("Number of exons/circRNA")
#Plot for Plasma-set
ggplot(count_plasma, aes(x = freq)) +
  geom_bar(fill = "orange") +
  theme_minimal() +
  ylab("freq") +
  xlab("# exons") +
  ggtitle("Number of exons/circRNA")
#Prepare combined dataset
full_counts <- bind_rows(count,
                         count_plasma,
                         .id = "set")
full_counts$set[full_counts$set == 1] <- "NB_set"
full_counts$set[full_counts$set == 2] <- "Plasma_set"
#Plot combined set
ggplot(full_counts, aes(x = freq, colour = set)) +
  geom_bar() +
  theme_minimal() +
  ylab("freq") +
  xlab("# exons") +
  ggtitle("Number of exons/circRNA")

###LENGTHS ---------------------------------------------------------------------
#Load in circRNA length data
circ_lengths <- read.delim("LOCATION OF CIRCRNA_LENGTHS", header=FALSE)
circ_lengths$V2 <- as.numeric(circ_lengths$V2)
#Plot for NB-set
ggplot(circ_lengths, aes(x = V2)) +
  geom_density() +
  theme_minimal() +
  xlab("# nucleotides") +
  ggtitle("Lengths of circRNAs")
#Plot for NB-set (x-axis reduced)
ggplot(circ_lengths, aes(x = V2)) +
  geom_density() +
  theme_minimal() +
  xlab("# nucleotides") +
  ggtitle("Lengths of circRNAs") +
  coord_cartesian(xlim=c(0,5000))
#Isolate circRNA lengths for plasma circRNAs
circ_plasma_lengths <- circ_lengths %>%
  filter(V1 %in% circ_plasma)
#Plot for Plasma-set
ggplot(circ_plasma_lengths, aes(x = V2)) +
  geom_density() +
  theme_minimal() +
  xlab("# nucleotides") +
  ggtitle("Lengths of circRNAs")
#Plot for Plasma-set (x-axis reduced)
ggplot(circ_plasma_lengths, aes(x = V2)) +
  geom_density() +
  theme_minimal() +
  xlab("# nucleotides") +
  ggtitle("Lengths of circRNAs") +
  coord_cartesian(xlim=c(0,5000))
#Prepare combined dataset
lengths <- bind_rows(circ_lengths,
                     circ_plasma_lengths,
                     .id = "set")
lengths[lengths == 1] <- "NB_set"
lengths[lengths == 2] <- "Plasma_set"
#Plot combined set
ggplot(lengths, aes(x = V2, color = set)) +
  geom_density() +
  theme_minimal() +
  xlab("# nucleotides") +
  ggtitle("Lengths of circRNAs")

###END -------------------------------------------------------------------------