#================================= PIPELINE RNA-seq KD vs Control =================================#

#---------------------------------PARAMETERS------------------------------------

n_top_rows <- 2000   #numbers of transcripts to be used

#---------------------------------LIBRARIES-------------------------------------

library(tidymodels)
library(ranger)
library(vip)
library(matrixStats)

#--------------------------------DATA LOADING-----------------------------------

dataset <- "Datasets/GSE266566_full.csv"
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

#run PCA to rank genes for feature selection
pca_result <- prcomp(t(log_top), center = TRUE, scale. = TRUE)

#------------------ML: CLASSIFICATION USING RANDOM FOREST-----------------------
set.seed(222) #reproducible results

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

#add knockdown target col from metadata
rf_df$knockdown <- factor(metadata$knockdown, levels = c("Control", "KD"))
rf_df$cell_line <- factor(metadata$cell_line)

#Knockdown recognition model

#split our df to train the model (75% train / 25% test)
splitkd  <- initial_split(rf_df, prop = 0.75, strata = knockdown)
trainkd  <- training(split)
testkd   <- testing(split)

#running Random Forest
rf_spec_kd <- rand_forest(mode = "classification", trees = 700) %>%
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

#Cell_line recognition model

#split our df to train the model (75% train / 25% test)
split_cl  <- initial_split(rf_df, prop = 0.75, strata = knockdown)
train_cl  <- training(split_cl)
test_cl   <- testing(split_cl)

#running Random Forest
rf_spec_cl <- rand_forest(mode = "classification", trees = 700) %>%
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
