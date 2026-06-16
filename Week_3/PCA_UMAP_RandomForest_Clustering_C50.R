#================================= PIPELINE RNA-seq KD vs Control =================================#

#---------------------------------PARAMETERS------------------------------------

n_top_rows <- 2000   #numbers of transcripts to be used

#---------------------------------LIBRARIES-------------------------------------

library(tidymodels)
#loading tidymodels loads the following libraries:
# --> library(ggplot2)
# --> library(rsample) for initial_split()
# --> library(parsnip) for rand_forest(), set_engine(), extract_fit_parsnip()
# --> library(recipes) for recipe(), step_zv()
# --> library(workflows) for workflow(), add_model(), add_recipe()
library(dplyr)
library(rsample) #for initial_split() --> Random Forest, C5.0
library(ranger)
library(matrixStats) #for rowVars()
library(vip) #for vip() in Random Forest Model

#--------------------------------DATA LOADING-----------------------------------

#import dataset from a .csv file
dataset <- "GSE266566_COUNTS.rsem_human.transcripts.csv"
df_BRCA <- read.csv2(dataset)

#remove variables that aren't going to be reused
rm(dataset)

#--------------------------------CLEANING---------------------------------------

#rename columns from column 1,2... to actual names
colnames(df_BRCA) <- df_BRCA[1,] 
df_BRCA <- df_BRCA[-1,]

#delete useless ensemble IDs
df_BRCA$`Ensembl.114.Transcript.ID` <- NULL

#convert column from chr to numeric to use them
df_BRCA[,3:54] <- lapply(df_BRCA[,3:54], function(x) as.numeric(as.character(x)))

#-----------------------------------LABELS--------------------------------------                   

# 52 samples (use the names of the column as sample names)
sample_names <- colnames(df_BRCA)[3:54] 

metadata <- data.frame(
  sample = sample_names,
  cell_line = sub("counts\\.([^.]+)\\..*", "\\1", sample_names),
  condition = sub("counts\\.[^.]+\\.([^.]+\\.[^.]+)\\.RNA.*", "\\1", sample_names),
  stringsAsFactors = FALSE
)

metadata$knockdown <- ifelse(grepl("NRP1", metadata$condition), "KD", "Control")
metadata$tech      <- ifelse(grepl("^sh", metadata$condition), "shRNA", "siRNA")

#remove variables that aren't going to be reused
rm(sample_names)

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

#remove variables that aren't going to be reused
rm(expr_mat, keep, lib_sizes, cpm_mat, vars, log_cpm, top_idx)



#--------------------------------FULL DATA PCA----------------------------------
#load librairies
library(FactoMineR)
library(factoextra)

#run PCA prompt
pca_result <- prcomp(t(log_top), center = TRUE, scale. = TRUE)

#PCA results on the 8 first dimensions
pca_data  <- as.data.frame(pca_result$x[, 1:8])
pca_data  <- cbind(pca_data, metadata) 
pca_var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

#PCA plot on PC1 & PC2
ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA — NRP1 knockdown vs Control",
       x = paste0("PC1 (", pca_var_explained[1], "%)"),
       y = paste0("PC2 (", pca_var_explained[2], "%)")) +
  theme_minimal()

#PCA plot on PC3 & PC4
ggplot(pca_data, aes(x = PC3, y = PC4, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA — NRP1 knockdown vs Control",
       x = paste0("PC3 (", pca_var_explained[3], "%)"),
       y = paste0("PC4 (", pca_var_explained[4], "%)")) +
  theme_minimal()

#prompt for accessing eigenvalues
# --> eigenvalue = percentage of explained variance by a dimension
pca_eig_val <- get_eigenvalue(pca_result)

#bar plot of the eigenvalues 
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))
#at dimension 6, we have 81.6% of the variance explained == keep 6 dimensions for machine learning 

#prompt for accessing results
pca_var <- get_pca_var(pca_result)

#most contributing descriptors in PC1 (first dim)
fviz_contrib(pca_result, choice = "var", axes = 1, top = 40)

#extract loadings from PC1 and PC2 and put them in dataframe
pca_loadings <- data.frame(
  Variable = rownames(pca_result$rotation),
  PC1 = pca_result$rotation[,1],
  PC2 = pca_result$rotation[,2]
)

#sort drivers from PC1 and PC2
pc1_loadings_sorted <- pca_loadings[order(abs(pca_loadings$PC1), decreasing = TRUE), ]
pc2_loadings_sorted <- pca_loadings[order(abs(pca_loadings$PC2), decreasing = TRUE), ]

#remove variables that are not useful for other methods
rm(pca_data, pca_var_explained, pca_eig_val, pca_var, pca_loadings, pc1_loadings_sorted, pc2_loadings_sorted)

#------------------------------------UMAP---------------------------------------
library(umap)

#UMAP
umap_res <- umap(t(log_top))
umap_df  <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df  <- cbind(umap_df, metadata)

#UMAP plot
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "UMAP — NRP1 knockdown vs Control") +
  theme_minimal()

#remove variables that aren't going to be reused
rm(umap_res,umap_df)

#-----------------------------PCA LOADINGS EXTRACT------------------------------
# extract loadings 
loadings <- pca_result$rotation  # genes × PCs

# choose Principal Components to consider, we take the first 6 --> see eigenvalues
pcs_to_use <- 1:6

# compute importance score per gene
gene_scores <- rowSums(abs(loadings[, pcs_to_use]))

# select top genes
top_gene_n <- min(500, length(gene_scores))
top_genes <- names(sort(gene_scores, decreasing = TRUE))[1:top_gene_n]

#create a dataframe with only the top genes from log_top
#transpose df to have our genes as cols
rf_df <- as.data.frame(t(log_top[top_genes, ]))

rm(loadings, pcs_to_use, gene_scores, top_gene_n)

#------------------ML: CLASSIFICATION USING RANDOM FOREST-----------------------
#set seed for reproducible results
set.seed(222)

#add knockdown target col from metadata
rf_df$knockdown <- factor(metadata$knockdown, levels = c("Control", "KD"))
rf_df$cell_line <- factor(metadata$cell_line)


#-----------------Knockdown recognition model

#split our df to train the model (75% train / 25% test)
split_kd  <- initial_split(rf_df, prop = 0.75, strata = knockdown)
train_kd  <- training(split_kd)
test_kd   <- testing(split_kd)

#running Random Forest using ranger engine and permutation importance
# --> permutation importance  in the ranger model for random forests measures how much the model's prediction error
#     increases when the values of a feature are shuffled
#     == helps id which features are most important for the model's predictions

rf_spec_kd <- rand_forest(mode = "classification", trees = 700, mtry = floor(sqrt(ncol(train_kd)-1))) %>%
  set_engine("ranger", importance = "permutation")

rf_recipe_kd <- recipe(knockdown ~ ., data = train_kd) %>%
  step_zv(all_predictors())

#fitting the model --> reduces the risk of overfitting, averaging the predictions of multiple trees
rf_workflow_kd <- workflow() %>%
  add_model(rf_spec_kd) %>%
  add_recipe(rf_recipe_kd)

rf_fit_kd <- fit(rf_workflow_kd, data = train_kd)

#---Plotting Random Forest results---

#predictions
preds_kd <- predict(rf_fit_kd, test_kd) %>%
  bind_cols(test_kd)

#confusion matrix --> plotted as a heatmap
conf_mat(preds_kd, truth = knockdown, estimate = .pred_class) %>% autoplot(type = "heatmap")

#model performance
metrics(preds_kd, truth = knockdown, estimate = .pred_class)

#variable importance from the fitted ranger model
rf_model_kd <- extract_fit_parsnip(rf_fit_kd)$fit
rf_importance_kd <- vip(rf_model_kd, num_features = 20)

#plot of the variable importance
rf_importance_kd

#remove variables that aren't going to be reused
rm(split_kd, train_kd, test_kd, rf_spec_kd, rf_recipe_kd, rf_workflow_kd, rf_fit_kd, preds_kd, rf_model_kd, rf_importance_kd)


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

#fitting the model
rf_workflow_cl <- workflow() %>%
  add_model(rf_spec_cl) %>%
  add_recipe(rf_recipe_cl)

rf_fit_cl <- fit(rf_workflow_cl, data = train_cl)

#---Plotting Random Forest results---

#predictions
preds_cl <- predict(rf_fit_cl, test_cl) %>%
  bind_cols(test_cl)

#confusion matrix --> plotted as a heatmap
conf_mat(preds_cl, truth = cell_line, estimate = .pred_class) %>% autoplot(type = "heatmap")

#model performance
metrics(preds_cl, truth = cell_line, estimate = .pred_class)

#variable importance from the fitted ranger model
rf_model_cl <- extract_fit_parsnip(rf_fit_cl)$fit
rf_importance_cl <- vip(rf_model_cl, num_features = 20)

#plot of the variable importance
rf_importance_cl

#remove variables that aren't going to be reused
rm(split_cl, train_cl, test_cl, rf_spec_cl, rf_recipe_cl, rf_workflow_cl, rf_fit_cl, preds_cl, rf_model_cl, rf_importance_cl)



#--------------------------HIERARCHICAL CLUSTERING------------------------------
#load libraries
library(factoextra)
library(cluster) #for agnes()

#--------cell lines using genes

#drop cell_line and knockdown columns --> agnes needs a numeric matrix
df_clust_cl <- rf_df[, !names(rf_df) %in% c("knockdown", "cell_line")]

#scale dataframe since the genes have different variance ranges
df_clust_cl <- as.data.frame(scale(df_clust_cl))

#define linkage methods
methods <- c( "average", "single", "complete", "ward")
names(methods) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
agglo_coef <- function(x) {
  agnes(df_clust_cl, method = x)$ac
}

#calculate agglomerative coefficient for each clustering linkage method
#the closer this value is to 1, the stronger the clusters
sapply(methods, agglo_coef)

#perform hierarchical clustering using chosen method (here ward is better)
clust_cl <- agnes(df_clust_cl, method = "ward")

#produce dendrogram
pltree(clust_cl, cex = 0.6, hang = -1, main = "Dendrogram") 

#remove variables that aren't going to be reused
rm(df_clust_cl, methods, agglo_coef, clust_cl)


#--------genes using cell lines

#scale dataframe since the genes have different variance ranges
df_clust_kd <- as.data.frame(scale(log_top[top_genes, ]))

#drop cell_line and knockdown columns --> agnes needs a numeric matrix
df_clust_kd <- df_clust_kd[, !names(df_clust_kd) %in% c("knockdown", "cell_line")]

#define linkage methods
methods <- c( "average", "single", "complete", "ward")
names(methods) <- c( "average", "single", "complete", "ward")

#function to compute agglomerative coefficient
agglo_coef <- function(x) {
  agnes(df_clust_kd, method = x)$ac
}

#calculate agglomerative coefficient for each clustering linkage method
#the closer this value is to 1, the stronger the clusters
sapply(methods, agglo_coef)

#perform hierarchical clustering using chosen method (here ward is better)
clust_kd <- agnes(df_clust_kd, method = "ward")

#produce dendrogram
pltree(clust_kd, cex = 0.6, hang = -1, main = "Dendrogram") 

#remove variables that aren't going to be reused
rm(df_clust_kd, methods, agglo_coef, clust_kd)


#------------------------------------C5.0---------------------------------------
#load libraries
library(C50) #for C5.0(); as.party.C5.0()

#add cell line (what we cant to classify)
rf_df$cell_line <- factor(metadata$cell_line)

#split our df to train the model (75% train / 25% test)
split <- initial_split(rf_df, prop = 0.75, strata = cell_line)
train <- training(split)
test <- testing(split)

#create C 5.0 model
model <- C5.0(cell_line ~ ., data = train)

#
#Afficher un résumé textuel de l'arbre pour vérifier les noms des colonnes
summary(model)
party_model <- as.party.C5.0(model)

#decision tree graph
plot(model)

#remove variables that aren't going to be reused
rm(split, train, test, model, party_model)
