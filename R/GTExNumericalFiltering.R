library(tidyverse)
library(ComplexHeatmap)
library(Cairo)
library(magick)
library(circlize)
library(factoextra)
library(NbClust)
library(tidyHeatmap)
library(viridis)
library(matrixStats)
library(RColorBrewer)

#Set seed for later references (747)
set.seed(747)

###RESOURCES -------------------------------------------------------------------
tissue_expression <- read.delim(file = "LOCATION OF GTEX MEDIAN EXPRESSION",
                                skip=2)
circ_ens_set <- read.delim("LOCATION OF CIRCRNA_ENS_IDS FROM CONSENSUS SUBSET")
circ_subset <- unique(circ_ens_set$gene_id)
NB_expression <- read.delim(file = "LOCATION OF STAR EXPRESSION COUNTS",
                            stringsAsFactors=FALSE)
NB_sample_names <- scan("LOCATION OF STAR EXPRESSION COUNTS NAMES",
                        what = "",
                        quiet = TRUE)
sorted_subset <- read.csv("LOCATION OF FULL CIRCRNAS FROM CONSENSUS SUBSET",
                          sep="",
                          stringsAsFactors=FALSE)
NB_restricted_set_full <- read.delim("LOCATION OF FULL CONSENSUS SUBSET")


###LOAD IN NB SUBSET -----------------------------------------------------------
#Median TPMs
NBMedianTpm <- function(NB_expression,
                        subset_genes){
  #Subset interested circRNA host genes
  NB_expression_subset <- subset(NB_expression,
                                 Geneid %in% subset_genes)
  NB_subset <- NB_expression_subset[,-1:-6]
  rownames(NB_subset) <- NB_expression_subset$Geneid
  NB_matrix <- data.matrix(NB_subset)
  NB_subset$Neuroblastoma <- rowMedians(NB_matrix)
  NB_median_tpm <- dplyr::select(NB_subset,
                                 Neuroblastoma)
  return(NB_median_tpm)
}
NB_subset <- NBMedianTpm(NB_expression,
                         circ_subset)
NB_expression_subset <- subset(NB_expression,
                               Geneid %in% circ_subset)
NB_subset$Name <- NB_expression_subset$Geneid
#Individual TPMs
NB38Tpm <- function(NB_expression,
                    subset_genes){
  #Subset interested circRNA host genes
  NB_expression_subset <- subset(NB_expression,
                                 Geneid %in% subset_genes)
  NB_subset <- NB_expression_subset[,-1:-6]
  rownames(NB_subset) <- NB_expression_subset$Geneid
  colnames(NB_subset) <- NB_sample_names
  return(NB_subset)
}
NB_38 <- NB38Tpm(NB_expression,
                 circ_subset)
NB_38$Name <- NB_expression_subset$Geneid

###DATA ------------------------------------------------------------------------
#Split Ensembl-stable ID and transcript version number
tissue_name_split <- separate(tissue_expression,
                              col = Name,
                              into = c("Name","Version"),
                              sep = "[.]")
#Subset interested circRNA host genes and removing PAR-transcripts
tissue_expression_subset <- subset(tissue_name_split,
                                   Name %in% circ_subset)
tissue_expression_subset <- tissue_expression_subset[!grepl("PAR",
                                                            tissue_expression_subset$Version),]
tissue_expression_subset <- left_join(tissue_expression_subset,
                                      NB_38,
                                      by = "Name")
tissue_subset <- tissue_expression_subset[,-1:-3]
tissue_subset$Testis <- NULL
tissue_subset$PMABM000CFV <- NULL
tissue_subset$PMABM000EIN <- NULL
tissue_subset$PMABM000GNS <- NULL
rownames(tissue_subset) <- tissue_expression_subset$Name
#Generation of NBMedian
tissue_nonlog_NB <- tissue_subset[,54:88]
tissue_subset$NBMedian <- rowMedians(data.matrix(tissue_nonlog_NB))
tissue_subset <- tissue_subset[,-54:-88]
# tissue_nonlog_NB <- tissue_subset[,54:88]
# tissue_subset$NBMedian <- rowMaxs(data.matrix(tissue_nonlog_NB))
# tissue_subset <- tissue_subset[,-54:-88]
#Numerical filtering
tissue_nonlog <- filter(tissue_subset,
                        NBMedian > 4)
tissue_matrix_nonlog <- data.matrix(tissue_nonlog[,1:53])
tissue_nonlog$GTExMean <- rowMeans(tissue_matrix_nonlog)
tissue_nonlog <- filter(tissue_nonlog,
                               GTExMean < 2)
tissue_matrix_nonlog <- data.matrix(tissue_nonlog[,1:53])
tissue_nonlog$GTExMax <- rowMaxs(tissue_matrix_nonlog)
tissue_nonlog <- filter(tissue_nonlog,
                        GTExMax < NBMedian)
tissue_nonlog$Ratio <- tissue_nonlog$NBMedian/tissue_nonlog$GTExMax
tissue_nonlog <- filter(tissue_nonlog,
                        Whole.Blood < 1)
tissue_meta <- tissue_nonlog[,54:57]
tissue_nonlog <- tissue_nonlog[,-54:-57]
#Replace NBMedian by 35NBs
tissue_nonlog$Name <- rownames(tissue_nonlog)
tissue_nonlog_NB$Name <- rownames(tissue_nonlog_NB)
tissue_nonlog_full <- left_join(tissue_nonlog, tissue_nonlog_NB, by = "Name")
rownames(tissue_nonlog_full) <- tissue_nonlog_full$Name
tissue_nonlog_full <- tissue_nonlog_full[,-54]
#Summary
summary <- summary(tissue_nonlog_full)
#log10 + 1 to ensure correct human visualisation
tissue_matrix <- log10(data.matrix(tissue_nonlog_full) + 1)
tissue_matrix_nonlog <- data.matrix(tissue_nonlog_full)

###DENSITYTPM ------------------------------------------------------------------
#Log10+1 TPM
TPM_density <- as_tibble(c(tissue_matrix))
ggplot(TPM_density, aes(x = value)) +
  geom_density(kernel = "gaussian")
median(TPM_density$value)

###CONVERSION GENE NAMES -------------------------------------------------------
last_ht <- as.data.frame(tissue_matrix)
NB_set_gene_names <- NB_restricted_set_full[,2:3]
NB_set_gene_names <- unique(NB_set_gene_names)
last_ht$gene_id <- rownames(last_ht)
last_ht <- left_join(last_ht,
                     NB_set_gene_names,
                     by = "gene_id")
rownames(last_ht) <- last_ht$gene_name
last_ht <- as.matrix(last_ht[,1:88])

###HEATMAP ---------------------------------------------------------------------
#Generate heatmap
ht <- Heatmap(last_ht,
              use_raster = T,
              row_dend_reorder = T,
              show_column_names = F,
              show_row_names = T,
              show_column_dend = F,
              show_row_dend = T,
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              row_split = 2,
              row_title = "CircRNA Host-Genes",
              column_split = 2,
              column_title = "Neuroblastoma (TPMs) vs GTEx (Median TPMs)",
              column_title_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 6),
              row_names_gp = gpar(fontsize = 8),
              col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
              heatmap_legend_param = list(title = "TPM (log10 + 1)",
              direction = "horizontal"))
#Draw heatmap to fix clusters
ht_drawn <- draw(ht,
                 heatmap_legend_side = "bottom")

###OUTPUT ----------------------------------------------------------------------
#Save heatmaps
save_pdf(ht,
         "OUTPUT LOCATION",
         units = "in",
         width = 8.3,
         height = 11.7)

save_pdf(ht,
         "OUTPUT LOCATION",
         units = "in",
         height = 8.3,
         width = 11.7)
         # width = 8.3,
         # height = 11.7)
numerical_subset <- tissue_nonlog_full[unlist(row_order(ht_drawn)),]
numerical_subset_names <- rownames(numerical_subset)
#Write tables to disk
write.table(numerical_subset_names,
            file = "OUTPUT LOCATION",
            quote = F,
            col.names = F,
            row.names = F)
write.table(numerical_subset,
            file = "OUTPUT LOCATION",
            quote = F,
            col.names = T,
            row.names = T)

###END -------------------------------------------------------------------------