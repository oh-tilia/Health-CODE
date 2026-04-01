#---------------------------------LOADING DATA----------------------------------

#Loading libraries
library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(umap)

#Dataset we will be using 
dataset <- "GSE266566_COUNTS.rsem_human.transcripts.csv"

#Loading a dataset from a csv file into a dataframe called df_BRCA
df_BRCA <- read.csv2(dataset)

#Remove variables that are no longer needed
rm(dataset)

#------------------------------DATA CLEANING------------------------------------

#Rename columns 
colnames(df_BRCA) <- df_BRCA[1,]
df_BRCA <- df_BRCA[-1,]

#Delete Ensembl ID
df_BRCA$`Ensembl.114.Transcript.ID` <- NULL

#Convert chr type cols to num type for PCA
df_BRCA[,3:54] <- lapply(df_BRCA[,3:54], function(x) as.numeric(as.character(x)))

#Checking the type of the data and the dataframe structure
str(df_BRCA)

#---------------------------------PCA-------------------------------------------

#Select only necessary variables (omit name variables)
X <- df_BRCA[,3:54]

#Log2 transformation of counts (recommended for RNA-seq)
X_log <- log2(X + 1)

#Remove genes (rows) with zero variance
gene_var <- apply(X_log, 1, var)

#Keep top 2000 most variable genes
top_genes <- order(gene_var, decreasing = TRUE)[1:2000]
X_filtered <- X_log[top_genes, ]

#Transpose matrix (samples = rows)
X_t <- t(X_filtered)

#Standardize the data
scaled_X <- scale(X_t)

#Compute PCA
pca_result <- prcomp(scaled_X)

#Extract principal components for each observation
pca_data <- as.data.frame(pca_result$x[,1:2])
pca_data$Sample <- rownames(pca_data)

#Groups
pca_data$CellLine <- ifelse(grepl("HS578T", pca_data$Sample), "HS578T",
                            ifelse(grepl("MDAMB231", pca_data$Sample), "MDAMB231","SUM159"))

pca_data$Treatment <- ifelse(grepl("NT", pca_data$Sample), "Control","NRP1 knockdown")

var_explained <- (pca_result$sdev^2)/sum(pca_result$sdev^2)

#----------------------------------PCA PLOT-------------------------------------

ggplot(pca_data, aes(x = PC1, y = PC2, color = CellLine, shape = Treatment)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Sample), size = 3) +
  labs(
    title = "PCA of RNA-seq Data",
    x = paste0("PC1 (", round(var_explained[1]*100,2), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2]*100,2), "% variance)")
  ) +
  theme_minimal()

#------------------------HEATMAP-----------------------------------------------

annotation_col <- data.frame(CellLine = pca_data$CellLine,
                             Treatment = pca_data$Treatment)
rownames(annotation_col) <- pca_data$Sample

pheatmap(X_filtered, 
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         main = "Top 2000 variable genes Heatmap")