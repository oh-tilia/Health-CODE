#================================= PIPELINE RNA-seq KD vs Control =================================#

#------------------------------PARAMÈTRES------------------------------------
n_top_rows <- 20000   # Nombre de transcrits les plus variables à utiliser

#------------------------------LIBRARIES--------------------------------------
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

#------------------------------CHARGEMENT DATA---------------------------------
dataset <- "GSE266566_COUNTS.rsem_human.transcripts.csv"
df_BRCA <- read.csv2(dataset)
rm(dataset)

#------------------------------NETTOYAGE---------------------------------------
colnames(df_BRCA) <- df_BRCA[1,]
df_BRCA <- df_BRCA[-1,]
df_BRCA$`Ensembl.114.Transcript.ID` <- NULL
df_BRCA[,3:54] <- lapply(df_BRCA[,3:54], function(x) as.numeric(as.character(x)))

#------------------------------LABELS-----------------------------------------
sample_names <- colnames(df_BRCA)[3:54]  # 52 échantillons
metadata <- data.frame(
  sample = sample_names,
  cell_line = sub("counts\\.([^.]+)\\..*", "\\1", sample_names),
  condition = sub("counts\\.[^.]+\\.([^.]+\\.[^.]+)\\.RNA.*", "\\1", sample_names),
  stringsAsFactors = FALSE
)
metadata$knockdown <- ifelse(grepl("NRP1", metadata$condition), "KD", "Control")
metadata$tech      <- ifelse(grepl("^sh", metadata$condition), "shRNA", "siRNA")

#------------------------------NORMALISATION-----------------------------------
expr_mat <- as.matrix(df_BRCA[, 3:54])
rownames(expr_mat) <- df_BRCA$`Ensembl.114.Transcript.Name`

# Filtre : transcrits exprimés dans au moins 3 échantillons
keep <- rowSums(expr_mat > 1) >= 3
expr_mat <- expr_mat[keep, ]
cat("Transcrits retenus après filtre:", nrow(expr_mat), "\n")

# Normalisation log2(CPM + 1)
lib_sizes <- colSums(expr_mat)
cpm_mat   <- sweep(expr_mat, 2, lib_sizes, "/") * 1e6
log_cpm   <- log2(cpm_mat + 1)

# Ajustement automatique du paramètre n_top_rows
if (n_top_rows > nrow(log_cpm)) {
  warning("n_top_rows > nombre de transcrits filtrés. Utilisation de ", nrow(log_cpm), " transcrits à la place.")
  n_top_rows <- nrow(log_cpm)
}

#------------------------------PCA & UMAP--------------------------------------
# Sélection des n_top_rows transcrits les plus variables
vars      <- rowVars(log_cpm)
top_idx   <- order(vars, decreasing = TRUE)[1:n_top_rows]
log_top   <- log_cpm[top_idx, ]

# PCA
pca_res <- prcomp(t(log_top), scale. = TRUE)
pca_df  <- as.data.frame(pca_res$x[, 1:5])
pca_df  <- cbind(pca_df, metadata)
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "PCA — NRP1 knockdown vs Control",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal()

# UMAP
umap_res <- umap(t(log_top))
umap_df  <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df  <- cbind(umap_df, metadata)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cell_line, shape = knockdown)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(title = "UMAP — NRP1 knockdown vs Control") +
  theme_minimal()

#------------------------------ML: CLASSIFICATION (k-fold CV pour tuning)-----------------------------
set.seed(42)

# Split global (75% train / 25% test)
split  <- initial_split(ml_df %>% select(-cell_line), prop = 0.75, strata = label)
train  <- training(split)
test   <- testing(split)

# Recipe avec normalisation et up-sampling pour équilibrer les classes
rec <- recipe(label ~ ., data = train) |>
  step_zv(all_predictors()) |>
  step_normalize(all_predictors()) |>
  step_upsample(label)

# Random Forest avec tuning
rf_spec <- rand_forest(trees = 500, mtry = tune(), min_n = tune()) |>
  set_engine("ranger", importance = "impurity") |>
  set_mode("classification")

wf <- workflow() |>
  add_recipe(rec) |>
  add_model(rf_spec)

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