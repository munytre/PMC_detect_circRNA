library(tidyverse)
library(DESeq2)
library(tximportData)
library(tximport)
library(biomaRt)
library(GenomicFeatures)
library(AnnotationHub)
library(ensembldb)
library(apeglm)
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
library(EnhancedVolcano)

#Set seed for later references (747)
set.seed(747)

###DATA ------------------------------------------------------------------------
#Load in path of NBs
NB_outfiles <- list.files(path = "LOCATION OF SALMON OUTPUT",
                          pattern = ".sf", 
                          full.names = T, 
                          recursive = T)
#Jobs file is a file containing sample names of the sequencing data
NB_sample_names <- scan("LOCATION OF JOBS FILE (SALMON)",
                        what = "",
                        quiet=TRUE)
names(NB_outfiles) <- NB_sample_names
#Load in path of GTEx
Salmon_outfiles <-  list.files(path = "LOCATION OF SALMON OUTPUT",
                               pattern = ".sf", 
                               full.names = T, 
                               recursive = T)
#Jobs file is a file containing sample names of the sequencing data
sample_names <- scan("LOCATION OF JOBS FILE (SALMON)",
                     what = "",
                     quiet=TRUE)
names(Salmon_outfiles) <- sample_names
#Merge paths
complete_outfiles <- c(NB_outfiles, Salmon_outfiles)
names(complete_outfiles)

###IMPORT QUANT VALUES ---------------------------------------------------------
#Import Salmon quant files, only transcript-level values
txi.complete.tx <- tximport(complete_outfiles,
                            type = "salmon",
                            txOut = T)

###GENERATE TX2GENE(v102) ------------------------------------------------------
#Make table to map transcript IDs to gene IDs for gene-level DE analysis
ah <- AnnotationHub()
query(ah, "EnsDb.HSapiens")
edb <- ah[["AH89180"]]
txs <- as_tibble(transcripts(edb,
                             return.type = "DataFrame"))
gxs <- as_tibble(genes(edb, return.type = "DataFrame"))
gxs_names <- gxs[,1:2]
txs_extended <- left_join(txs,
                          gxs_names,
                          by = "gene_id")
tx2gene <- dplyr::select(txs_extended,
                         tx_id_version,
                         gene_name)
###SUMMARIZETOGENE -------------------------------------------------------------
#Summarize transcript-level expression to gene-level expression (NBs)
txi.complete.gene <- summarizeToGene(txi.complete.tx,
                               tx2gene = tx2gene)
#Subset estimated genecounts as data frame
counts.complete <- as.data.frame(txi.complete.gene$counts)

###EXPORT GENE ESTIMATIONS -----------------------------------------------------
#Export for gene estimations
# write.table(counts.complete,
#             file = "OUTPUT LOCATION",
#             sep = "\t",
#             quote = F,
#             row.names = T,
#             col.names = T)

###DESEQ2 (METADATA) -----------------------------------------------------------
#Fix metadata for DESeq2
complete_metadata <- as_tibble(names(complete_outfiles))
complete_metadata <- complete_metadata %>%
  mutate(source = value)
complete_metadata[1:35,2] <- "Tumor"
complete_metadata[36:594,2] <- "GTEx"
colnames(complete_metadata) <- c("sample_name",
                                 "source")
GTEx_metadata <- read.delim("LOCATION OF GTEX METADATA")
GTEx_sample_tissue <- as.data.frame(GTEx_metadata$tissue_type_detail)
GTEx_sample_tissue$specimen_id <- GTEx_metadata$specimen_id
colnames(GTEx_sample_tissue) <- c("tissue_type",
                                  "sample_name")
complete_metadata <- left_join(complete_metadata,
                               GTEx_sample_tissue,
                               by = "sample_name") %>%
  as.data.frame()
complete_metadata[is.na(complete_metadata)] <- "Neuroblastoma"
rownames(complete_metadata) <- complete_metadata$sample_name
complete_metadata <- complete_metadata[,-1]

###EXPORT METADATA -------------------------------------------------------------
#Export for gene counts metadata
# write.table(complete_metadata,
#             file = "OUTPUT LOCATION",
#             sep = "\t",
#             quote = F,
#             row.names = T,
#             col.names = T)

###DESEQ2 ----------------------------------------------------------------------
#Import txi.complete.gene in DESeq2
dds <- DESeqDataSetFromTximport(txi = txi.complete.gene,
                                colData = complete_metadata,
                                design = ~ source)
#Perform DESeq2
dds <- DESeq(dds)
#Isolate the counts generated in DESeq2 (res)
resultsNames(dds)
res <- results(dds,
               contrast = c("source", "Tumor", "GTEx"))
res05 <- results(dds,
                 alpha = 0.05,
                 contrast = c("source", "Tumor", "GTEx"))
#Order res by significance
resOrdered <- res[order(res$pvalue),]
#Subset the gene names that are significant and have a log2FC > 2 or -2
res_de_genes <- rownames(subset(resOrdered,
                                (resOrdered$padj < 0.05 & resOrdered$log2FoldChange > 1)))
#vst transformation (normalization) for the DESeq values
vsd <- vst(dds, blind=FALSE)

###DESEQ2 PCA ------------------------------------------------------------------
#DESeq2 PCA plots
plotPCA(vsd, intgroup=c("source"))
plotPCA(vsd, intgroup=c("source", "tissue_type"))

###DESEQ2 FILTER (NUMERICAL) ---------------------------------------------------
#Load in names of genes found in the numerical filtering (TPM)
numerical_filtering_genes <- as.data.frame(scan("LOCATION OF NAMES OF GENES FROM NUMERICAL FILTERING (TPM)",
                                                what = "",
                                                quiet=TRUE))
#Convert full set to gene_names
colnames(numerical_filtering_genes) <- "gene_id"
numerical_filtering_genes_set <- left_join(numerical_filtering_genes,
                      gxs,
                      by = "gene_id")
numerical_filtering_genes <- c(numerical_filtering_genes_set$gene_name)
#Subsetting genes found in TPM heatmap for validation
# high <- numerical_filtering_genes[37:56]
# numerical_filtering_genes_high <- res_de_genes[res_de_genes %in% high]
numerical_filtering_genes_DE <- res_de_genes[res_de_genes %in% numerical_filtering_genes]
ht_counts <- assay(vsd, normalized = T)[numerical_filtering_genes_DE,]
#Save DE genes
# write.table(numerical_filtering_genes_DE,
#             file = "OUTPUT LOCATION",
#             sep = "\t",
#             quote = F,
#             row.names = F,
#             col.names = F)

###CONVERSION GENE NAMES -------------------------------------------------------
NB_restricted_set_full <- read.delim("LOCATION OF FULL CONSENSUS SUBSET")
last_ht <- as.data.frame(ht_counts)
last_ht <- t(scale(t(last_ht)))
last_ht <- as.matrix(last_ht)

###COMPLEXHEATMAPS (NUMERICAL) -------------------------------------------------
#Set heatmap parameters
ht <- Heatmap(last_ht,
              use_raster = T,
              row_dend_reorder = T,
              show_row_names = T,
              show_row_dend = F,
              show_column_dend = F,
              show_column_names = F,
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              # row_split = 2,
              column_split = 4,
              column_gap = unit(3, "mm"),
              row_title = "CircRNA Host-Genes",
              column_title = "Neuroblastoma vs GTEx (DESeq2)",
              column_title_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 6),
              col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(256),
              heatmap_legend_param = list(title = "Scaled expression"),
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#BA313D","#31BAAE","#757575","#31BAAE"),
                                                                            col = c("#BA313D","#31BAAE","#757575","#31BAAE")))))
#Draw the heatmap
ht_drawn <- draw(ht)
#Saving the heatmap
save_pdf(ht,
         "OUTPUT LOCATION",
         units = "in",
         width = 11.7,
         height = 8.3)

###CONVERSION GENE NAMES (circs) -----------------------------------------------
last_ht_circ <- as.data.frame(ht_counts)
last_ht_circ$gene <- rownames(last_ht_circ)
last_ht_circ <- filter(last_ht_circ,
                       gene %in% c("CCDC7",
                                   "RANBP17", 
                                   "MBD5", 
                                   "ZNF124",
                                   "SPATA5"))
last_ht_circ <- last_ht_circ[,-595]
last_ht_circ <- t(scale(t(last_ht_circ)))
last_ht_circ <- as.matrix(last_ht_circ)

###COMPLEXHEATMAPS (circs) -----------------------------------------------------
#Set heatmap parameters
ht_circ <- Heatmap(last_ht_circ,
              use_raster = T,
              row_dend_reorder = T,
              show_row_names = T,
              show_row_dend = F,
              show_column_dend = F,
              show_column_names = F,
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              # row_split = 2,
              column_split = 2,
              column_gap = unit(3, "mm"),
              row_title = "CircRNA Host-Genes",
              column_title = "Neuroblastoma vs GTEx (DESeq2)",
              column_title_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 6),
              col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(256),
              heatmap_legend_param = list(title = "Scaled expression"),
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#BA313D","#31BAAE"),
                                                                            col = c("#BA313D","#31BAAE")))))
#Draw the heatmap
ht_drawn_circ <- draw(ht_circ)
#Saving the heatmap
save_pdf(ht_circ,
         "OUTPUT LOCATION",
         units = "in",
         width = 11.7,
         height = 8.3)

###COMPLEXHEATMAPS (NUMERICAL NON-SCALED) --------------------------------------
#Set heatmap parameters
ht <- Heatmap(ht_counts,
              use_raster = T,
              row_dend_reorder = T,
              show_row_names = T,
              show_row_dend = F,
              show_column_dend = F,
              show_column_names = F,
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              # row_split = 2,
              column_split = 8,
              column_gap = unit(3, "mm"),
              row_title = "CircRNA Host-Genes",
              column_title = "Neuroblastoma vs GTEx (DESeq2)",
              column_title_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 6),
              col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(256),
              heatmap_legend_param = list(title = "Relative expression"),
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#BA313D",
                                                                                     "#BA313D",
                                                                                     "#31BAAE",
                                                                                     "#31BAAE",
                                                                                     "#31BAAE",
                                                                                     "#31BAAE",
                                                                                     "#31BAAE",
                                                                                     "#31BAAE"),
                                                                            col = c("#BA313D",
                                                                                    "#BA313D",
                                                                                    "#31BAAE",
                                                                                    "#31BAAE",
                                                                                    "#31BAAE",
                                                                                    "#31BAAE",
                                                                                    "#31BAAE",
                                                                                    "#31BAAE")))))
#Draw the heatmap
ht_drawn <- draw(ht)
#Saving the heatmap
save_pdf(ht,
         "OUTPUT LOCATION",
         units = "in",
         width = 11.7,
         height = 8.3)

### TISSUE TYPE TEST -----------------------------------------------------------
test3 <- last_ht_circ[,column_order(ht_drawn_circ)[[2]]]
test3_1 <- as.data.frame(colnames(test3))
View(test3_1)
View(GTEx_metadata)
View(GTEx_sample_tissue)
colnames(test3_1) <- "sample_name"
test3_2 <- left_join(test3_1, GTEx_sample_tissue, by = "sample_name")

###ENHANCED VOLCANO (NUMERICAL) ------------------------------------------------
#Load in host genes with plasma circRNAs
DE_genes_circs <- read.delim("LOAD IN DE GENES WITH CIRCRNA IN PLASMA")
#Subset of genes with found circRNAs
found_circs <- as.data.frame(unique(DE_genes_circs$gene_id))
colnames(found_circs) <- "gene_id"
found_circs <- left_join(found_circs,
                         gxs,
                         by = "gene_id")
found_circs <- c(found_circs$gene_name)
#Setting colors for gene of interest
keyvals.colour <- ifelse(
  rownames(res) %in% found_circs, 'blue',
  ifelse(
    rownames(res) %in% numerical_filtering_genes_DE, 'red',
    'grey'))
keyvals.colour[is.na(keyvals.colour)] <- 'grey'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'Others'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'Neuroblastoma restricted genes (subset)'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Neuroblastoma restricted genes (subset) with plasma circRNAs'
#Setting shape for gene of interest
keyvals.shape <- ifelse(
  rownames(res) %in% found_circs, 16,
  ifelse(
    rownames(res) %in% numerical_filtering_genes_DE, 16,
    46))
keyvals.shape[is.na(keyvals.shape)] <- 46
names(keyvals.shape)[keyvals.shape == 46] <- 'Others'
names(keyvals.shape)[keyvals.shape == 16] <- 'Neuroblastoma restricted genes (subset)'
names(keyvals.shape)[keyvals.shape == 16] <- 'Neuroblastoma restricted genes (subset) with plasma circRNAs'
#Making the volcano plot
EnhancedVolcano(res,
                title = "Differential expressed genes (Neuroblastoma vs GTEx)",
                lab = rownames(res),
                colCustom = keyvals.colour,
                pointSize = c(ifelse(rownames(res) %in% numerical_filtering_genes_DE,
                                     3,
                                     0.5)),
                xlim = c(-10,10),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                pCutoffCol = 'padj',
                drawConnectors = T,
                directionConnectors = 'both',
                shapeCustom = keyvals.shape,
                selectLab = found_circs,
                legendPosition = 'top')

###CHECK ENSG PRESENCE ---------------------------------------------------------
#Check if ENSG is present in dds object
ENSG_check <- as_tibble(setdiff(numerical_filtering_genes,
                                numerical_filtering_genes_DE))
ENSG_check <- ENSG_check %>%
  mutate(in_dds = value %in% rownames(dds))

nonDE_genes <- ENSG_check %>%
  dplyr::filter(in_dds == TRUE) %>%
  pull(value)

###PLOTS COUNTS (NUMERICAL) ----------------------------------------------------
#Plot the normalized counts for genes of interest (DE)
pdf("OUTPUT LOCATION")
for(i in numerical_filtering_genes_DE){
  plotC <- plotCounts(dds,
                      gene = i,
                      intgroup="source",
                      transform = F,
                      normalized = T)
  print(plotC)
}
dev.off()
#Plot the normalized counts for genes of interest (NON_DE)
pdf("OUTPUT LOCATION")
for(i in nonDE_genes){
  plotC <- plotCounts(dds,
                      gene = i,
                      intgroup="source",
                      transform = F,
                      normalized = T)
  print(plotC)
}
dev.off()

###CONSENSUS NB SUBSET ---------------------------------------------------------
#Load in NB-restricted host-genes
NB_restricted_set_full <- read.delim("LOCATION OF FULL CONSENSUS SUBSET")
NB_set_gene_names <- NB_restricted_set_full[,2:3]
NB_set_gene_names <- unique(NB_set_gene_names)
circ_subset <- c(NB_set_gene_names$gene_id)
#Make data frame from DESeq2 res
DE_res_df <- as_tibble(res, rownames = NA)
#DE_res_df_plasma <- DE_res_df[found_circs,]
DE_res_df_plasma <- dplyr::filter(DE_res_df,
                                  rownames(DE_res_df) %in% found_circs)
#Subset the gene names that are significant and have a log2FC > 2 or -2
res_de_genes <- rownames(subset(resOrdered,
                                (resOrdered$padj < 0.05 & resOrdered$log2FoldChange > 4)))
#Filter the NB-restricted host-genes from the DE-subset
NB_restricted_de_genes <- res_de_genes[res_de_genes %in% circ_subset]
ht_counts <- assay(vsd, normalized = T)[NB_restricted_de_genes,]
#Save DE genes (NB-RESTRICTED SUBSET)
write.table(NB_restricted_de_genes,
            file = "OUTPUT LOCATION",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)

###ENHANCED VOLCANO (CONSENSUS NB SUBSET) --------------------------------------
#Making the volcano plot
EnhancedVolcano(res,
                title = "Differential expressed genes (Neuroblastoma vs GTEx)",
                lab = rownames(res),
                pointSize = c(ifelse(rownames(res) %in% NB_restricted_de_genes,
                                     3,
                                     0.5)),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-100,
                pCutoffCol = 'padj',
                selectLab = found_circs,
                legendPosition = 'top')

###COMPLEXHEATMAPS (CONSENSUS NB SUBSET) ---------------------------------------
#Set heatmap parameters
ht <- Heatmap(ht_counts,
              use_raster = T,
              row_dend_reorder = T,
              show_row_names = T,
              show_row_dend = F,
              show_column_dend = F,
              show_column_names = T,
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              row_split = 6,
              column_split = 8,
              row_title = "Host genes",
              column_title = "GTEx vs NBs",
              column_title_gp = gpar(fontsize = 15),
              column_names_gp = gpar(fontsize = 6),
              col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
              heatmap_legend_param = list(title = "Relative expression"))
#Draw the heatmap
ht_drawn <- draw(ht)

###FULL INFO ON GENES ----------------------------------------------------------
full_set <- left_join(DE_genes_circs,
                      gxs,
                      by = "gene_id")
DE_res_df_plasma$gene_name <- rownames(DE_res_df_plasma)
full_set <- left_join(full_set,
                      DE_res_df_plasma,
                      by = "gene_name")
report_table <- full_set[,c(1:5,18,22)]
report_table <- report_table[order(report_table$padj),]
#Save full set
write.table(full_set,
            file = "OUTPUT LOCATION",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
#Save report table
write.table(report_table,
            file = "OUTPUT LOCATION",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

### TEST NB MARKERS ------------------------------------------------------------
res_test[res_test %in% c('CHRNA3',
                         'DBH',
                         'FMO3',
                         'GAP43',
                         'GUSB',
                         'PHOX2B',
                         'POSTN',
                         'PRRX1',
                         'TH')]

###END -------------------------------------------------------------------------