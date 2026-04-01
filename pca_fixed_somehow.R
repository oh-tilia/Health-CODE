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

#Select only necessary variables (omit name variables).
rownames(df_BRCA) <- df_BRCA[,1]
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
pca_data <- as.data.frame(pca_result$x[,1:4])
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


#PCA 2
ggplot(pca_data, aes(x = PC3, y = PC4, color = CellLine, shape = Treatment)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Sample), size = 3) +
  labs(
    title = "PCA of RNA-seq Data",
    x = paste0("PC1 (", round(var_explained[1]*100,2), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2]*100,2), "% variance)")
  ) +
  theme_minimal()

#------------------------------------HEATMAP------------------------------------

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


#-----------------------------PCA DRIVERS---------------------------------------

# Charger les packages nécessaires
library(ggplot2)
library(tidyverse)
library("FactoMineR")
library("factoextra")

# run PCA prompt
# don't center all variables cause some shouldn't be centered (placette/Placet and id)
pca_result <- prcomp(X_filtered, center = FALSE, scale. = FALSE)

eig_val <- get_eigenvalue(pca_result)

#bar plot of the eigenvalue (percentage of explained variance by that dimension)
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))

# prompt for accessing results
var <- get_pca_var(pca_result)

# most contributing descriptors in PCA1 (first dim)
var_pca1_sort <- sort(var$contrib[,1], decreasing = TRUE)

# visualize PCA variables
fviz_pca_var(pca_result, col.var = "pink")

# visualize PCA variables V1
fviz_pca_var(pca_result, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

# visualize PCA variables V2
fviz_pca_var(pca_result, col.var = "contrib") +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 1) +
  theme_minimal()

# visualize PCA individuals
fviz_pca_ind(pca_result, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
#Montre l'importance des individus dans la PCA

fviz_pca_biplot(pca_result, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

# extract loadings from PC1 and PC2
loadingsPC1 <- pca_result$rotation[,1]
loadingsPC2 <- pca_result$rotation[,2]

# put them in dataframe
loadingsPCA <- data.frame(
  Variable = rownames(pca_result$rotation),
  PC1 = pca_result$rotation[,1],
  PC2 = pca_result$rotation[,2]
)

# sort drivers
loadingsPC1_sorted <- loadingsPCA[order(abs(loadingsPCA$PC1), decreasing = TRUE), ]
loadingsPC2_sorted <- loadingsPCA[order(abs(loadingsPCA$PC2), decreasing = TRUE), ]


rm(pca_data,pca_result,X)


#-------------------------PCA 2 GPT-----------------------
#---------------------------------LOADING DATA----------------------------------

library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(factoextra)

dataset <- "GSE266566_COUNTS.rsem_human.transcripts.csv"
df_BRCA <- read.csv2(dataset)

rm(dataset)

#------------------------------DATA CLEANING------------------------------------

# Rename columns
colnames(df_BRCA) <- df_BRCA[1,]
df_BRCA <- df_BRCA[-1,]

# ⚠️ IMPORTANT : garder les noms de gènes
# 👉 adapte "gene_name" si besoin (mets le bon nom de colonne)
rownames(df_BRCA) <- df_BRCA$gene_name

# Remove unnecessary columns
df_BRCA$`Ensembl.114.Transcript.ID` <- NULL
df_BRCA$gene_name <- NULL

# Convert to numeric
df_BRCA[,] <- lapply(df_BRCA[,], function(x) as.numeric(as.character(x)))

str(df_BRCA)

#---------------------------------PCA-------------------------------------------

# Matrix
X <- as.matrix(df_BRCA)

# Log transform
X_log <- log2(X + 1)

# Variance filtering
gene_var <- apply(X_log, 1, var)

top_genes <- order(gene_var, decreasing = TRUE)[1:2000]
X_filtered <- X_log[top_genes, ]

# Keep gene names
rownames(X_filtered) <- rownames(X_log)[top_genes]

# Transpose → samples = rows
X_t <- t(X_filtered)

# Scale
scaled_X <- scale(X_t)

# PCA
pca_result <- prcomp(scaled_X, center = TRUE, scale. = TRUE)

# PCA data
pca_data <- as.data.frame(pca_result$x[,1:4])
pca_data$Sample <- rownames(pca_data)

# Groups
pca_data$CellLine <- ifelse(grepl("HS578T", pca_data$Sample), "HS578T",
                            ifelse(grepl("MDAMB231", pca_data$Sample), "MDAMB231","SUM159"))

pca_data$Treatment <- ifelse(grepl("NT", pca_data$Sample), "Control","NRP1 knockdown")

# Variance explained
var_explained <- (pca_result$sdev^2)/sum(pca_result$sdev^2)

#----------------------------------PCA PLOT-------------------------------------

ggplot(pca_data, aes(x = PC1, y = PC2, color = CellLine, shape = Treatment)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Sample), size = 3) +
  labs(
    title = "PCA of RNA-seq Data",
    x = paste0("PC1 (", round(var_explained[1]*100,2), "%)"),
    y = paste0("PC2 (", round(var_explained[2]*100,2), "%)")
  ) +
  theme_minimal()

# PCA 2
ggplot(pca_data, aes(x = PC3, y = PC4, color = CellLine, shape = Treatment)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Sample), size = 3) +
  labs(
    title = "PCA of RNA-seq Data",
    x = paste0("PC3 (", round(var_explained[3]*100,2), "%)"),
    y = paste0("PC4 (", round(var_explained[4]*100,2), "%)")
  ) +
  theme_minimal()

#------------------------------------HEATMAP------------------------------------

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

#-----------------------------PCA DRIVERS---------------------------------------

# Variance explained plot
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))

# Loadings = contribution des gènes
loadings <- pca_result$rotation

loadingsPCA <- data.frame(
  Gene = rownames(loadings),
  PC1 = loadings[,1],
  PC2 = loadings[,2]
)

# Top drivers
top_PC1 <- loadingsPCA %>%
  arrange(desc(abs(PC1))) %>%
  slice(1:20)

top_PC2 <- loadingsPCA %>%
  arrange(desc(abs(PC2))) %>%
  slice(1:20)

print(top_PC1)
print(top_PC2)

# PCA variables (gènes)
fviz_pca_var(pca_result,
             col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# PCA individus (samples)
fviz_pca_ind(pca_result,
             col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

# Biplot
fviz_pca_biplot(pca_result,
                repel = TRUE,
                col.var = "#2E9FDF",
                col.ind = "#696969")

#------------------------------------CLEAN--------------------------------------

#rm(X, X_log, X_filtered, X_t, scaled_X)
