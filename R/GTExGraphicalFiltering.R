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

#Set seed for later references (747)
set.seed(747)

###RESOURCES -------------------------------------------------------------------
tissue_expression <- read.delim(file = "LOCATION OF GTEX MEDIAN EXPRESSION",
                                skip=2)
circ_subset <- scan("LOCATION OF HOSTGENES_ENS_IDS FROM CONSENSUS SUBSET",
                    what = "",
                    quiet=TRUE)
NB_expression <- read.delim(file = "LOCATION OF STAR EXPRESSION COUNTS",
                            stringsAsFactors=FALSE)
NB_sample_names <- scan("LOCATION OF STAR EXPRESSION COUNTS NAMES",
                        what = "",
                        quiet = TRUE)
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
  NB_median_tpm <- select(NB_subset,
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
#NB_subset or NB_38
tissue_expression_subset <- left_join(tissue_expression_subset,
                                      NB_subset,
                                      by = "Name")
tissue_subset <- tissue_expression_subset[,-1:-3]
tissue_subset$Testis <- NULL
tissue_subset$PMABM000CFV <- NULL
tissue_subset$PMABM000EIN <- NULL
tissue_subset$PMABM000GNS <- NULL
rownames(tissue_subset) <- tissue_expression_subset$Name
tissue_matrix_nonlog <- data.matrix(tissue_subset)
summary <- summary(tissue_matrix_nonlog)
#log10 + 1 to ensure correct human visualisation
tissue_matrix <- log10(data.matrix(tissue_subset) + 1)

###DENSITYTPM ------------------------------------------------------------------
#Non-log TPM
TPM_density <- as_tibble(c(tissue_matrix_nonlog))
ggplot(TPM_density, aes(x = value)) +
  geom_density(kernel = "gaussian")
median(TPM_density$value)
#Log10+1 TPM
TPM_density <- as_tibble(c(tissue_matrix))
ggplot(TPM_density, aes(x = value)) +
  geom_density(kernel = "gaussian")
median(TPM_density$value)

###HEATMAP ---------------------------------------------------------------------
#Generate heatmap
ht <- Heatmap(tissue_matrix,
              use_raster = T,
              row_dend_reorder = T,
              show_row_names = F,
              show_row_dend = F,
              clustering_method_rows = "average",
              row_split = 12,
              row_title = "CircRNA Host-Genes",
              column_split = 12,
              column_title = "Neuroblastoma (Median TPMs) vs GTEx (Median TPMs)",
              column_title_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 6),
              heatmap_legend_param = list(title = "TPM (log10 + 1)",
                                          at = c(0,2)),
              col = inferno(256))
#Draw heatmap to fix clusters
ht_drawn <- draw(ht)
#BLT = Blood, lymphocytes and testis
save_pdf(ht,
         "OUTPUT LOCATION",
         units = "in",
         width = 8.3,
         height = 11.7)
#Isolate genes from the generated clusters
ht_cluster7 <- t(t(row.names(tissue_subset[row_order(ht_drawn)[[7]],])))
ht_cluster8 <- t(t(row.names(tissue_subset[row_order(ht_drawn)[[8]],])))
ht_cluster9 <- t(t(row.names(tissue_subset[row_order(ht_drawn)[[9]],])))
ht_cluster10 <- t(t(row.names(tissue_subset[row_order(ht_drawn)[[10]],])))
ht_cluster11 <- t(t(row.names(tissue_subset[row_order(ht_drawn)[[11]],])))
ht_cluster12 <- t(t(row.names(tissue_subset[row_order(ht_drawn)[[12]],])))
ht_cluster_total <- c(ht_cluster7,
                      ht_cluster8,
                      ht_cluster9,
                      ht_cluster10,
                      ht_cluster11,
                      ht_cluster12)

###SUBSETHT --------------------------------------------------------------------
#Subset the genes (1st heatmap) of interest from the main data set
tissue_ht_subset <- subset(tissue_expression_subset,
                                   Name %in% ht_cluster_total)
tissue_ht <- tissue_ht_subset[,-1:-3]
tissue_ht$Testis <- NULL
rownames(tissue_ht) <- tissue_ht_subset$Name
ht_matrix_nonlog <- data.matrix(tissue_ht)
#log10 + 1 to ensure correct human visualisation
ht_matrix <- log10(data.matrix(tissue_ht) + 1)

#Log10+1 TPM
TPM_density_ht <- as_tibble(c(ht_matrix))
ggplot(TPM_density_ht, aes(x = value)) +
  geom_density(kernel = "gaussian")
median(TPM_density_ht$value)

#Generate silhouette plots to determine numbers of k in k-means clustering
#Columnwise
fviz_nbclust(ht_matrix,
             kmeans,
             method = "silhouette") +
  labs(subtitle = "Silhouette method (Original subset) - Transformed")

#Generate heatmap
subset_ht <- Heatmap(ht_matrix,
              use_raster = T,
              row_dend_reorder = T,
              show_row_names = F,
              show_row_dend = F,
              clustering_method_rows = "average",
              row_km = 4,
              # row_split = 6,
              row_title = "CircRNA Host-Genes",
              column_split = 6,
              column_title = "Neuroblastoma (Median TPMs) vs GTEx (Median TPMs)",
              column_title_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 6),
              heatmap_legend_param = list(title = "TPM (log10 + 1)",
                                          at = c(0,2)),
              col = inferno(256))
#Draw heatmap to fix clusters
subset_ht_drawn <- draw(subset_ht)
save_pdf(subset_ht,
         "OUTPUT LOCATION",
         units = "in",
         width = 8.3,
         height = 11.7)
#Subset the genes (2st heatmap) of interest
km_matrix <- t(t(tissue_ht[row_order(subset_ht_drawn)[[3]],]))

###SUBSETKMEANS ----------------------------------------------------------------
#Log10+1 TPM
TPM_density_km <- as_tibble(c(km_matrix))
ggplot(TPM_density_km, aes(x = value)) +
  geom_density(kernel = "gaussian")
median(TPM_density_km$value)

#Generate silhouette plots to determine numbers of k in k-means clustering
#Columnwise
fviz_nbclust(km_matrix,
             kmeans,
             method = "silhouette") +
  labs(subtitle = "Silhouette method (Original subset) - Transformed")

#Generate heatmap
subset_km <- Heatmap(km_matrix,
                     use_raster = T,
                     row_dend_reorder = T,
                     show_row_names = F,
                     show_row_dend = F,
                     clustering_method_rows = "average",
                     row_km = 6,
                     row_title = "CircRNA Host-Genes",
                     column_split = 4,
                     column_title = "Neuroblastoma (Median TPMs) vs GTEx (Median TPMs)",
                     column_title_gp = gpar(fontsize = 15),
                     column_names_gp = gpar(fontsize = 6),
                     heatmap_legend_param = list(title = "TPM",
                                                 at = c(0,2)),
                     col = inferno(256))
#Draw heatmap to fix clusters
subset_km_drawn <- draw(subset_km)
save_pdf(subset_km,
         "OUTPUT LOCATION",
         units = "in",
         width = 8.3,
         height = 11.7)
last_ht <- as.data.frame(t(t(km_matrix[row_order(subset_km_drawn)[[3]],])))
#Assign genenames for last heatmap based on genes (3rd heatmap) of interest
NB_set_gene_names <- unique(NB_restricted_set_full[,2:3])
last_ht$gene_id <- rownames(last_ht)
last_ht <- left_join(last_ht,
                     NB_set_gene_names,
                     by = "gene_id")
rownames(last_ht) <- last_ht$gene_name
last_ht <- as.matrix(last_ht[,1:54])

###LASTHT ----------------------------------------------------------------------
#Log10+1 TPM
TPM_density_km <- as_tibble(c(last_ht))
ggplot(TPM_density_km, aes(x = value)) +
  geom_density(kernel = "gaussian")
median(TPM_density_km$value)

#Generate silhouette plots to determine numbers of k in k-means clustering
#Columnwise
fviz_nbclust(last_ht,
             kmeans,
             method = "silhouette") +
  labs(subtitle = "Silhouette method (Original subset) - Transformed")

#Generate heatmap
last_ht_ht <- Heatmap(last_ht,
                     use_raster = T,
                     row_dend_reorder = T,
                     show_row_names = T,
                     show_row_dend = F,
                     clustering_method_rows = "average",
                     # row_km = 2,
                     row_split = 6,
                     row_title = "CircRNA Host-Genes",
                     column_split = 2,
                     column_title = "Neuroblastoma (Median TPMs) vs GTEx (Median TPMs)",
                     column_title_gp = gpar(fontsize = 15),
                     column_names_gp = gpar(fontsize = 6),
                     heatmap_legend_param = list(title = "TPM",
                                                 at = c(0,50)),
                     col = inferno(256))
#Draw heatmap to fix clusters
last_ht_drawn <- draw(last_ht_ht)
save_pdf(last_ht_ht,
         "OUTPUT LOCATION",
         units = "in",
         width = 8.3,
         height = 11.7)
#Isolate tables containing genes of interest
sorted_subset <- last_ht[unlist(row_order(last_ht_drawn)),]
sorted_subset_names <- rownames(sorted_subset)
#Write tables to disk
write.table(sorted_subset_names,
            file = "OUTPUT LOCATION",
            quote = F,
            col.names = F,
            row.names = F)
write.table(sorted_subset,
            file = "OUTPUT LOCATION",
            quote = F,
            col.names = T,
            row.names = T)

###END -------------------------------------------------------------------------
