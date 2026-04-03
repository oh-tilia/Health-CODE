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

#automatically adjust n_top_row parameter 
#if number of transcripts post filtration is lower than 20000, use the number of transcripts instead as a variable
if (n_top_rows > nrow(log_cpm)) {
  warning("n_top_rows > nombre de transcrits filtrés. Utilisation de ", nrow(log_cpm), " transcrits à la place.")
  n_top_rows <- nrow(log_cpm)
}

#------------------------------------PCA----------------------------------------
#selection of the transcripts that explain the most the variance out of n_top_rows
vars      <- rowVars(log_cpm)
top_idx   <- order(vars, decreasing = TRUE)[1:n_top_rows]
log_top   <- log_cpm[top_idx, ]   #FINAL DF TO BE USED

#PCA
pca_res <- prcomp(t(log_top), scale. = TRUE)
pca_data  <- as.data.frame(pca_res$x[, 1:5]) #PCA results on the 5 first dimensions
pca_data  <- cbind(pca_data, metadata) 
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA — NRP1 knockdown vs Control",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal()

#---------------------------------PCA DRIVERS-----------------------------------

#load librairies
library(ggplot2)
library(tidyverse)
library("FactoMineR")
library("factoextra")

#run PCA prompt
pca_result <- prcomp(t(log_top), center = TRUE, scale. = TRUE)

eig_val <- get_eigenvalue(pca_result)

#bar plot of the eigenvalue (percentage of explained variance by that dimension)
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))
#at dimension 8, we have 73.0% of the variance explained == keep 8 dim for machine learning 

#prompt for accessing results
var <- get_pca_var(pca_result)

#most contributing descriptors in PCA1 (first dim)
var_pca1_sort <- sort(var$contrib[,1], decreasing = TRUE)

#visualize PCA variables V1 TO BE FIXED
fviz_pca_var(pca_result, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

#visualize PCA variables V2  TO BE FIXED
fviz_pca_var(pca_result, col.var = "contrib") +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 1) +
  theme_minimal()

#visualize PCA individuals (less pretty then first graph)
fviz_pca_ind(pca_result, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

#biplot individuals + variables
fviz_pca_biplot(pca_result, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "pink"  # Individuals color
)

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

#UMAP
umap_res <- umap(t(log_top))
umap_df  <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df  <- cbind(umap_df, metadata)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "UMAP — NRP1 knockdown vs Control") +
  theme_minimal()



#------------------ML: CLASSIFICATION USING RANDOM FOREST-----------------------
set.seed(222) #reproducible results

#load libraries
library(tidymodels)
library(ranger)
library(rfviz)
library(tidyverse)

# extract loadings 
loadings <- pca_res$rotation  # genes × PCs

# choose Principal Components to consider, we take the first 5
pcs_to_use <- 1:5

# compute importance score per gene
gene_scores <- rowSums(abs(loadings[, pcs_to_use]))

# select top genes
top_genes <- names(sort(gene_scores, decreasing = TRUE))[1:500]

#create a dataframe with only the top 200 genes from log_top
#transpose df to have our genes as cols
rf_df <- as.data.frame(t(log_top[top_genes, ]))

#add cell line col from metadata
rf_df$cell_line <- metadata$cell_line

#split our df to train the model (75% train / 25% test)
split  <- initial_split(rf_df, prop = 0.75, strata = cell_line)
train  <- training(split)
test   <- testing(split)

#running Random Forest avec tuning
rf_spec <- rand_forest(mode = "classification", engine = "ranger", trees = 700, mtry = tune(), min_n = tune())

#fit the model
rf_fit <- fit(rf_spec, cell_line ~ ., data = train)

#---Plotting Random Forest results---

library(yardstick)

# predictions
preds <- predict(rf_fit, test) %>%
  bind_cols(test)

# confusion matrix
conf_mat(preds, truth = cell_line, estimate = .pred_class)

conf_mat(preds, truth = cell_line, estimate = .pred_class) %>% autoplot(type = "heatmap")

#PLOT DOESNT WORK
#MODEL NEEDS TUNING









#shit AI code

# k-fold CV pour tuning (k = 5)
folds <- vfold_cv(train, v = 5, strata = label)

# Grid pour tuning
grid <- grid_regular(
  mtry(range = c(10, 100)),
  min_n(range = c(2, 8)),
  levels = 3
)

# Tuning
tune_res <- tune_grid(
  wf,
  resamples = folds,
  grid = grid,
  metrics = metric_set(accuracy, roc_auc, f_meas)
)

# Meilleurs paramètres
best_params <- select_best(tune_res, metric = "roc_auc")
final_wf    <- finalize_workflow(wf, best_params)

# Évaluation finale sur le test set
final_fit <- last_fit(final_wf, split)
collect_metrics(final_fit)
