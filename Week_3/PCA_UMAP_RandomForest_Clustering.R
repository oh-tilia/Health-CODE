#================================= PIPELINE RNA-seq KD vs Control =================================#

#---------------------------------PARAMETERS------------------------------------

n_top_rows <- 2000   #numbers of transcripts to be used

#---------------------------------LIBRARIES-------------------------------------

library(tidymodels)
library(ranger)
library(vip)
library(matrixStats)

#--------------------------------DATA LOADING-----------------------------------

dataset <- "GSE266566_COUNTS.rsem_human.transcripts.csv"
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

#automatically adjust n_top_row parameter 
#if number of transcripts post filtration is lower than 20000, use the number of transcripts instead as a variable
if (n_top_rows > nrow(log_cpm)) {
  warning("n_top_rows > nombre de transcrits filtrés. Utilisation de ", nrow(log_cpm), " transcrits à la place.")
  n_top_rows <- nrow(log_cpm)
}

#selection of the transcripts that explain the most the variance out of n_top_rows
vars      <- rowVars(log_cpm)
top_idx   <- order(vars, decreasing = TRUE)[1:n_top_rows]
log_top   <- log_cpm[top_idx, ]   #FINAL DF TO BE USED


#--------------------------------FULL DATA PCA----------------------------------
#load librairies
library(ggplot2)
library(tidyverse)
library("FactoMineR")
library("factoextra")

#run PCA prompt
pca_result <- prcomp(t(log_top), center = TRUE, scale. = TRUE)

#PCA results on the 8 first dimensions
pca_data  <- as.data.frame(pca_result$x[, 1:8])
pca_data  <- cbind(pca_data, metadata) 
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

#PCA plot on PC1 & PC2
ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA — NRP1 knockdown vs Control",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal()

#PCA plot on PC3 & PC4
ggplot(pca_data, aes(x = PC3, y = PC4, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA — NRP1 knockdown vs Control",
       x = paste0("PC3 (", var_explained[3], "%)"),
       y = paste0("PC4 (", var_explained[4], "%)")) +
  theme_minimal()

eig_val <- get_eigenvalue(pca_result)

#bar plot of the eigenvalue (percentage of explained variance by that dimension)
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))
#at dimension 8, we have 73.0% of the variance explained == keep 8 dim for machine learning 

#prompt for accessing results
var <- get_pca_var(pca_result)

#most contributing descriptors in PCA1 (first dim)
var_pca1_sort <- sort(var$contrib[,1], decreasing = TRUE)

#extract loadings from PC1 and PC2
loadingsPC1 <- pca_result$rotation[,1]
loadingsPC2 <- pca_result$rotation[,2]

#put them in dataframe
loadingsPCA <- data.frame(
  Variable = rownames(pca_result$rotation),
  PC1 = pca_result$rotation[,1],
  PC2 = pca_result$rotation[,2]
)

#sort drivers
loadingsPC1_sorted <- loadingsPCA[order(abs(loadingsPCA$PC1), decreasing = TRUE), ]
loadingsPC2_sorted <- loadingsPCA[order(abs(loadingsPCA$PC2), decreasing = TRUE), ]



#------------------------------------UMAP---------------------------------------
library(umap)

#selection of the transcripts that explain the most the variance out of n_top_rows
vars      <- rowVars(log_cpm)
top_idx   <- order(vars, decreasing = TRUE)[1:n_top_rows]
log_top   <- log_cpm[top_idx, ]   #FINAL DF TO BE USED

#UMAP
umap_res <- umap(t(log_top))
umap_df  <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df  <- cbind(umap_df, metadata)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "UMAP — NRP1 knockdown vs Control") +
  theme_minimal()



#-----------------------------PCA LOADINGS EXTRACT------------------------------
# extract loadings 
loadings <- pca_result$rotation  # genes × PCs

# choose Principal Components to consider, we take the first 5
pcs_to_use <- 1:5

# compute importance score per gene
gene_scores <- rowSums(abs(loadings[, pcs_to_use]))

# select top genes
top_gene_n <- min(500, length(gene_scores))
top_genes <- names(sort(gene_scores, decreasing = TRUE))[1:top_gene_n]

#create a dataframe with only the top genes from log_top
#transpose df to have our genes as cols
rf_df <- as.data.frame(t(log_top[top_genes, ]))

#------------------ML: CLASSIFICATION USING RANDOM FOREST-----------------------
set.seed(222) #reproducible results

#add knockdown target col from metadata
rf_df$knockdown <- factor(metadata$knockdown, levels = c("Control", "KD"))
rf_df$cell_line <- factor(metadata$cell_line)

#-----------------Knockdown recognition model

#split our df to train the model (75% train / 25% test)
splitkd  <- initial_split(rf_df, prop = 0.75, strata = knockdown)
trainkd  <- training(splitkd)
testkd   <- testing(splitkd)

#running Random Forest
rf_spec_kd <- rand_forest(mode = "classification", trees = 700, mtry = floor(sqrt(ncol(trainkd)-1))) %>%
  set_engine("ranger", importance = "permutation")

rf_recipe_kd <- recipe(knockdown ~ ., data = trainkd) %>%
  step_zv(all_predictors())

#fit the model
rf_workflow_kd <- workflow() %>%
  add_model(rf_spec_kd) %>%
  add_recipe(rf_recipe_kd)

rf_fit_kd <- fit(rf_workflow_kd, data = trainkd)

#---Plotting Random Forest results---

# predictions
preds <- predict(rf_fit_kd, testkd) %>%
  bind_cols(testkd)

# confusion matrix
conf_mat(preds, truth = knockdown, estimate = .pred_class)

conf_mat(preds, truth = knockdown, estimate = .pred_class) %>% autoplot(type = "heatmap")

# model performance
metrics(preds, truth = knockdown, estimate = .pred_class)

# variable importance from the fitted ranger model
rf_model_kd <- extract_fit_parsnip(rf_fit_kd)$fit
rf_importance_kd <- vip(rf_model_kd, num_features = 20)

rf_importance_kd

#------Cell_line recognition model

#split our df to train the model (75% train / 25% test)
split_cl  <- initial_split(rf_df, prop = 0.75, strata = cell_line)
train_cl  <- training(split_cl)
test_cl   <- testing(split_cl)

#running Random Forest
rf_spec_cl <- rand_forest(mode = "classification", trees = 700, mtry = floor(sqrt(ncol(train_cl)-1)) ) %>%
  set_engine("ranger", importance = "permutation")

rf_recipe_cl <- recipe(cell_line ~ ., data = train_cl) %>%
  step_zv(all_predictors())

#fit the model
rf_workflow_cl <- workflow() %>%
  add_model(rf_spec_cl) %>%
  add_recipe(rf_recipe_cl)

rf_fit_cl <- fit(rf_workflow_cl, data = train_cl)

#---Plotting Random Forest results---

# predictions
preds_cl <- predict(rf_fit_cl, test_cl) %>%
  bind_cols(test_cl)

# confusion matrix
conf_mat(preds_cl, truth = cell_line, estimate = .pred_class)

conf_mat(preds_cl, truth = cell_line, estimate = .pred_class) %>% autoplot(type = "heatmap")

# model performance
metrics(preds_cl, truth = cell_line, estimate = .pred_class)

# variable importance from the fitted ranger model
rf_model_cl <- extract_fit_parsnip(rf_fit_cl)$fit
rf_importance_cl <- vip(rf_model_cl, num_features = 20)

rf_importance_cl


#--------------------------HIERARCHICAL CLUSTERING------------------------------
library(factoextra)
library(cluster)

#drop cell_line and knockdown columns --> agnes needs a numeric matrix
df_clust <- rf_df[, !names(rf_df) %in% c("knockdown", "cell_line")]

#scale dataframe since the genes have different variance ranges
df_clust <- as.data.frame(scale(df_clust))

#define linkage methods
methods <- c( "average", "single", "complete", "ward")
names(methods) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
agglo_coef <- function(x) {
  agnes(df_clust, method = x)$ac
}

#calculate agglomerative coefficient for each clustering linkage method
#the closer this value is to 1, the stronger the clusters
sapply(methods, agglo_coef)

#perform hierarchical clustering using chosen method (here ward is better)
clust <- agnes(df_clust, method = "ward")

#produce dendrogram
pltree(clust, cex = 0.6, hang = -1, main = "Dendrogram") 







library(factoextra)
library(cluster)

#scale dataframe since the genes have different variance ranges
df_clust <- as.data.frame(scale(log_top[top_genes, ]))

#drop cell_line and knockdown columns --> agnes needs a numeric matrix
df_clust <- df_clust[, !names(df_clust) %in% c("knockdown", "cell_line")]

#define linkage methods
methods <- c( "average", "single", "complete", "ward")
names(methods) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
agglo_coef <- function(x) {
  agnes(df_clust, method = x)$ac
}

#calculate agglomerative coefficient for each clustering linkage method
#the closer this value is to 1, the stronger the clusters
sapply(methods, agglo_coef)

#perform hierarchical clustering using chosen method (here ward is better)
clust <- agnes(df_clust, method = "ward")

#produce dendrogram
pltree(clust, cex = 0.6, hang = -1, main = "Dendrogram") 














