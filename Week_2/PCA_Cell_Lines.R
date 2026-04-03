#================================= PIPELINE RNA-seq KD vs Control =================================#

#---------------------------------PARAMETERS------------------------------------

n_top_rows <- 20000   #numbers of transcripts to be used

#---------------------------------LIBRARIES-------------------------------------

library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(umap)
library(tidymodels)
library(ranger)
library(vip)
library(matrixStats)
library(themis)  # pour step_upsample()

#--------------------------------DATA LOADING-----------------------------------

dataset <- "GSE266566_full.csv"
df_BRCA <- read.csv2(dataset)
rm(dataset)

#--------------------------------CLEANING---------------------------------------

#rename columns from column 1,2... to actual names
colnames(df_BRCA) <- df_BRCA[1,] 
df_BRCA <- df_BRCA[-1,]

#delete useless ensemble IDs
df_BRCA$`Ensembl.114.Transcript.ID` <- NULL

#convert column from chr to numeric to use them
df_BRCA[,3:54] <- lapply(df_BRCA[,3:54], function(x) as.numeric(as.character(x)))

#--------------------------------LABELS-----------------------------------------                   

#52 samples (use the names of the column as sample names)
sample_names <- colnames(df_BRCA)[3:54] 

#creating metadata
metadata <- data.frame(
  sample = sample_names,
  cell_line = sub("counts\\.([^.]+)\\..*", "\\1", sample_names),
  condition = sub("counts\\.[^.]+\\.([^.]+\\.[^.]+)\\.RNA.*", "\\1", sample_names),
  stringsAsFactors = FALSE
)

metadata$knockdown <- ifelse(grepl("NRP1", metadata$condition), "KD", "Control")
metadata$tech      <- ifelse(grepl("^sh", metadata$condition), "shRNA", "siRNA")

#-------------------------------NORMALIZATION-----------------------------------

#select only necessary variables (omit name variables)
expr_mat <- as.matrix(df_BRCA[, 3:54])

#use the name of  the transcript as indexes                         
rownames(expr_mat) <- df_BRCA$`Ensembl.114.Transcript.Name`

#keep transcripts that are present in AT LEAST 3 samples
keep <- rowSums(expr_mat > 1) >= 3
expr_mat <- expr_mat[keep, ]
cat("Transcrits retenus après filtre:", nrow(expr_mat), "\n")

#normalization w/ log2(CPM + 1)
lib_sizes <- colSums(expr_mat)
cpm_mat   <- sweep(expr_mat, 2, lib_sizes, "/") * 1e6
log_cpm   <- log2(cpm_mat + 1)

log_HS578T <- (data.frame(log_cpm) %>% select(contains("HS578T")))
log_MDAMB231 <- (data.frame(log_cpm) %>% select(contains("MDAMB231")))
log_SUM159 <- (data.frame(log_cpm) %>% select(contains("SUM159")))




#-----------------------------------PCA_all-------------------------------------
#selection of the transcripts that explain the most the variance out of n_top_rows

#automatically adjust n_top_row parameter 
#if number of transcripts post filtration is lower than 20000, use the number of transcripts instead as a variable
if (n_top_rows > nrow(log_cpm)) {
  warning("n_top_rows > nombre de transcrits filtrés. Utilisation de ", nrow(log_cpm), " transcrits à la place.")
  n_top_rows <- nrow(log_cpm)
}

vars      <- rowVars(log_cpm)
top_idx   <- order(vars, decreasing = TRUE)[1:n_top_rows]
log_top   <- log_cpm[top_idx, ]   #FINAL DF TO BE USED

pca_res <- prcomp(t(log_top), scale. = TRUE)
pca_data  <- as.data.frame(pca_res$x[, 1:8]) #PCA results on the 5 first dimensions
pca_data  <- cbind(pca_data, metadata) 
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA — NRP1 knockdown vs Control",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal()

#-----------------------------------PCA_HS578T----------------------------------


if (n_top_rows > nrow(log_HS578T)) {
  warning("n_top_rows > nombre de transcrits filtrés. Utilisation de ", nrow(log_cpm), " transcrits à la place.")
  n_top_rows <- nrow(log_HS578T)
}


vars      <- rowVars(data.matrix(log_HS578T))
top_idx   <- order(vars, decreasing = TRUE)[1:n_top_rows]
log_top   <- log_HS578T[top_idx, ]   #FINAL DF TO BE USED

#PCA_HS578T
pca_res <- prcomp(t(log_top), scale. = TRUE)
pca_data  <- as.data.frame(pca_res$x[, 1:8]) #PCA results on the 5 first dimensions
pca_data  <- cbind(pca_data, metadata[is.element(metadata$cell_line, "HS578T"), ]) 
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA — NRP1 knockdown vs Control _ Cell line HS578T",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal()

#-----------------------------------PCA_MDAMB231--------------------------------


if (n_top_rows > nrow(log_MDAMB231)) {
  warning("n_top_rows > nombre de transcrits filtrés. Utilisation de ", nrow(log_cpm), " transcrits à la place.")
  n_top_rows <- nrow(log_MDAMB231)
}

vars      <- rowVars(data.matrix(log_MDAMB231))
top_idx   <- order(vars, decreasing = TRUE)[1:n_top_rows]
log_top   <- log_MDAMB231[top_idx, ]   #FINAL DF TO BE USED

#PCA_MDAMB231
pca_res <- prcomp(t(log_top), scale. = TRUE)
pca_data  <- as.data.frame(pca_res$x[, 1:5]) #PCA results on the 5 first dimensions
pca_data  <- cbind(pca_data, metadata[is.element(metadata$cell_line, "MDAMB231"), ]) 
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

ggplot(pca_data, aes(x = PC1, y = PC2, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85, colour = "#7CAE00") +
  labs(title = "PCA — NRP1 knockdown vs Control _ Cell line MDAMB231",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal()

#-----------------------------------PCA_SUM159----------------------------------

if (n_top_rows > nrow(log_SUM159)) {
  warning("n_top_rows > nombre de transcrits filtrés. Utilisation de ", nrow(log_cpm), " transcrits à la place.")
  n_top_rows <- nrow(log_SUM159)
}

vars      <- rowVars(data.matrix(log_SUM159))
top_idx   <- order(vars, decreasing = TRUE)[1:n_top_rows]
log_top   <- log_SUM159[top_idx, ]   #FINAL DF TO BE USED

#PCA_MDAMB231
pca_res <- prcomp(t(log_top), scale. = TRUE)
pca_data  <- as.data.frame(pca_res$x[, 1:5]) #PCA results on the 5 first dimensions
pca_data  <- cbind(pca_data,metadata[is.element(metadata$cell_line, "SUM159"), ]) 
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

ggplot(pca_data, aes(x = PC1, y = PC2, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85,color = "#00BFC4") +
  labs(title = "PCA — NRP1 knockdown vs Control _ Cell line SUM159",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal()
