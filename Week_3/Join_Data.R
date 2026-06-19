#==============================================================================#
#                        PIPELINE Cancerous vs Healthy                         #
#==============================================================================#


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

#---------------------------------PARAMETERS------------------------------------

n_top_rows <- 2000   #numbers of transcripts to be used

#import data 
df_BRCA <- read.csv2("GSE266566_COUNTS.csv")
df_HD <-read.csv2("GSE71862_MCF7_MCF10A_RSEM_expectedcounts.csv")

#--------------------------------CLEANING---------------------------------------

#rename columns from column 1,2... to actual names
if (colnames(df_BRCA)[1] == 'Column1'){
  colnames(df_BRCA) <- df_BRCA[1,] 
  df_BRCA <- df_BRCA[-1,]
}

#rename columns from column 1,2... to actual names
if (colnames(df_HD)[1] == 'Column1'){
  colnames(df_HD) <- df_HD[1,] 
  df_HD <- df_HD[-1,]
}
#remove unneeded columns from the dataset
df_BRCA <- df_BRCA[ ,c(-1,-3)]
df_HD$accession <- NULL

#convert column from chr to numeric in order to use them
df_BRCA[,2:53] <- lapply(df_BRCA[,2:53], function(x) as.numeric(as.character(x)))
df_HD[,2:7] <- lapply(df_HD[,2:7], function(x) as.numeric(as.character(x)))

#--------------------------TREATMENT AND COMBINATION-----------------------------

#regex expression to get only the gene name for the transcript name                      
df_BRCA$gene <- sub('([^-]+).*','\\1',df_BRCA$Ensembl.114.Transcript.Name)

#total gene expression of the gene across all sample
# --> if duplicate genes, only the one with the most expression is kept
df_BRCA$total <-  rowSums(df_BRCA[,c(2:53)])

#replace the column with transcript name by the column with gene name so thats its on the left of the dataset
#maybe not the optimal way to do it but it works fine
#slice_max bottleneck
df_BRCA_reduced <- slice_max(df_BRCA, order_by= total , by=gene)
df_BRCA_reduced$Ensembl.114.Transcript.Name <- df_BRCA_reduced$gene
df_BRCA_reduced$gene <-NULL
names(df_BRCA_reduced)[names(df_BRCA_reduced) == 'Ensembl.114.Transcript.Name'] <- 'gene'

#combine the two dataset used to get a bigger one with 3 healthy samples and 50+ cancerous samples                      
df_combined <- inner_join(df_BRCA_reduced,df_HD)
df_combined$total <- NULL

#remove variables that aren't going to be reused
rm(df_BRCA,df_BRCA_reduced,df_HD)

#-----------------------------------LABELS--------------------------------------                   

# 58 samples (use the names of the column as sample names)
sample_names <- colnames(df_combined)[2:59] 

#create a metadata dataframe to easily access the cell line and state of each sample
metadata <- data.frame(
  sample = sample_names,
  cell_line = gsub('counts\\.|\\.s.*|_.*',"",sample_names),
  stringsAsFactors = FALSE
)

#add state (cancerous or healthy) of cells to the metadata
metadata$state <- ifelse(grepl("MCF10A", metadata$cell_line), "Healthy", "Cancerous")

#remove variables that aren't going to be reused
rm(sample_names)

#-------------------------------NORMALIZATION-----------------------------------

#select only necessary variables (omit name variables)
expr_mat <- as.matrix(df_combined[, 2:59])

#use the name of  the transcript as indexes                         
rownames(expr_mat) <- df_combined$gene

#keep transcripts that are present in AT LEAST 3 samples
keep <- rowSums(expr_mat > 1) >= 3
expr_mat <- expr_mat[keep, ]
cat("Transcrits retenus après filtre:", nrow(expr_mat), "\n")

#normalization w/ log2(CPM + 1)
#Log2 is used when normalizing the expression of genes because it aids in calculating fold change, 
#which measures the up-regulated vs down-regulated genes between samples. Log2 measured data is also 
#closer to the biologically-detectable changes.

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
#load librairies for PCA
library(FactoMineR)
library(factoextra)

#run PCA prompt
pca_result <- prcomp(t(log_top), center = TRUE, scale. = TRUE)

#PCA results on the 8 first dimensions
pca_data  <- as.data.frame(pca_result$x[, 1:8])
pca_data  <- cbind(pca_data, metadata) 
pca_var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

#PCA plot on PC1 & PC2
ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_line, shape = state)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA - cancerous vs healthy",
       x = paste0("PC1 (", pca_var_explained[1], "%)"),
       y = paste0("PC2 (", pca_var_explained[2], "%)")) +
  theme_minimal()

#PCA plot on PC3 & PC4
ggplot(pca_data, aes(x = PC3, y = PC4, color = cell_line, shape = state)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA - cancerous vs healthy",
       x = paste0("PC3 (", pca_var_explained[3], "%)"),
       y = paste0("PC4 (", pca_var_explained[4], "%)")) +
  theme_minimal()

#prompt for accessing eigenvalues
# --> eigenvalue = percentage of explained variance by a dimension
pca_eig_val <- get_eigenvalue(pca_result)

#bar plot of the eigenvalues 
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))
#at dimension 6, we have 88.4% of the variance explained == keep 6 dimensions for machine learning 

#the following commented lines of code can be useful to visualize the data, but it is not mandatory 
#fviz_contrib(...) doesnt work but since it is not necessary for the rest of the code to work it was not fixed
                      
# #prompt for accessing results
# pca_var <- get_pca_var(pca_result)
# 
# #most contributing descriptors in PC1 (first dim)
# fviz_contrib(pca_result, choice = 'var', axes = 1, top = 40)
# 
# #extract loadings from PC1 and PC2 and put them in dataframe
# pca_loadings <- data.frame(
#   Variable = rownames(pca_result$rotation),
#   PC1 = pca_result$rotation[,1],
#   PC2 = pca_result$rotation[,2]
# )
# 
# #sort drivers from PC1 and PC2
# pc1_loadings_sorted <- pca_loadings[order(abs(pca_loadings$PC1), decreasing = TRUE), ]
# pc2_loadings_sorted <- pca_loadings[order(abs(pca_loadings$PC2), decreasing = TRUE), ]
# 

#remove variables that are not useful for other methods
rm(pca_data, pca_var_explained, pca_eig_val, pca_var, pca_loadings, pc1_loadings_sorted, pc2_loadings_sorted)


#-----------------------------PCA LOADINGS EXTRACT------------------------------
#extract loadings 
loadings <- pca_result$rotation  # genes × PCs

#choose Principal Components to consider, we take the first 6 --> see eigenvalues
pcs_to_use <- 1:6

#compute importance score per gene
gene_scores <- rowSums(abs(loadings[, pcs_to_use]))

# select top genes
top_gene_n <- min(500, length(gene_scores))
top_genes <- names(sort(gene_scores, decreasing = TRUE))[1:top_gene_n]

#create a dataframe with only the top genes from log_top
#transpose df to have our genes as cols
rf_df <- as.data.frame(t(log_top[top_genes, ]))

#remove unneeded variables
rm(loadings, pcs_to_use, gene_scores, top_gene_n)

#------------------ML: CLASSIFICATION USING RANDOM FOREST-----------------------
#set seed for reproducible results
set.seed(222)

#add state target col from metadata
rf_df$state <- factor(metadata$state)
rf_df$cell_line <- factor(metadata$cell_line)

#-----------------Knockdown recognition model

#split our df to train the model (75% train / 25% test)
split  <- initial_split(rf_df, prop = 0.75, strata = state)
train  <- training(split)
test   <- testing(split)

#removing columns with duplicate data in the train and test
train <- train[, -c(361, 362, 374)]
test  <- test[,  -c(361, 362, 374)]

#running Random Forest using ranger engine and permutation importance
# --> permutation importance  in the ranger model for random forests measures how much the model's prediction error
#     increases when the values of a feature are shuffled
#     == helps id which features are most important for the model's predictions

rf_spec <- rand_forest(mode = "classification", trees = 700, mtry = floor(sqrt(ncol(train)-1))) %>%
  set_engine("ranger", importance = "permutation")

rf_recipe <- recipe(state ~ ., data = train) %>%
  step_zv(all_predictors())

#fitting the model --> reduces the risk of overfitting, averaging the predictions of multiple trees
rf_workflow <- workflow() %>%
  add_model(rf_spec) %>%
  add_recipe(rf_recipe)

rf_fit <- fit(rf_workflow, data = train)

#prints the underlying ranger object with all call details
rf_fit %>% extract_fit_parsnip() %>% .$fit

#---Plotting Random Forest results---

#predictions
preds <- predict(rf_fit, test) %>%
  bind_cols(test)

#confusion matrix --> plotted as a heatmap
conf_mat(preds, truth = state, estimate = .pred_class) %>% autoplot(type = "heatmap")

#model performance
metrics(preds, truth = state, estimate = .pred_class)

#variable importance from the fitted ranger model
rf_model <- extract_fit_parsnip(rf_fit)$fit
rf_importance <- vip(rf_model, num_features = 20)

#plot of the variable importance
rf_importance

#remove variables that aren't going to be reused
rm(split, train, test, rf_spec, rf_recipe, rf_workflow, rf_fit, preds, rf_model, rf_importance)


#------------------------------------C5.0---------------------------------------
#load libraries for C5.0
library(C50) #for C5.0(); as.party.C5.0()

#add cell line (what we cant to classify)
rf_df$cell_line <- factor(metadata$cell_line)

#split our df to train the model (75% train / 25% test)
split <- initial_split(rf_df, prop = 0.75, strata = cell_line)
train <- training(split)
test <- testing(split)

#removing columns with duplicate data in the train and test
train <- train[, -c(361, 362, 374)]
test  <- test[,  -c(361, 362, 374)]

#create C 5.0 model
model <- C5.0(cell_line ~ ., data = train)

#print text summary of the tree --> verify genes used
summary(model)
party_model <- as.party.C5.0(model)

#decision tree graph
plot(model)

#remove variables that aren't going to be reused
rm(split, train, test, model, party_model)

#------------------------------------C5.0---------------------------------------
#load libraries for C5.0
library(C50) #for C5.0(); as.party.C5.0()

#set seed for reproducible results
set.seed(666)

#remove cell_line 
# --> if not removed the model uses the cell_line label to separate states and does not rely on the genes
rf_df$cell_line <- NULL

#add state (what we cant to classify)
rf_df$state <- factor(metadata$state)

#split our df to train the model (75% train / 25% test)
split <- initial_split(rf_df, prop = 0.75, strata = state)
train <- training(split)
test <- testing(split)

#removing columns with duplicate data in the train and test
train <- train[, -c(361, 362, 374)]
test  <- test[,  -c(361, 362, 374)]

#create C 5.0 model
model <- C5.0(state ~ ., data = train)

#print text summary of the tree --> verify genes used
summary(model)
party_model <- as.party.C5.0(model)

#decision tree graph
plot(model)

#remove variables that aren't going to be reused
rm(split, train, test, model, party_model)


