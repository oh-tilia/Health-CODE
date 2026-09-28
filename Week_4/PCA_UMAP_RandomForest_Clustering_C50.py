#================================= PIPELINE RNA-seq KD vs Control =================================#
# Python 3.12 translation of PCA_UMAP_RandomForest_Clustering_C50.R
#
# Requirements:
#   pip install numpy pandas scipy scikit-learn matplotlib seaborn umap-learn

#---------------------------------PARAMETERS------------------------------------

n_top_rows = 2000   # numbers of transcripts to be used
SEED = 222          # seed for reproducible results

#---------------------------------LIBRARIES-------------------------------------

import re
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier, export_text, plot_tree
from sklearn.model_selection import train_test_split
from sklearn.metrics import (accuracy_score, cohen_kappa_score,
                             ConfusionMatrixDisplay)
import umap

sns.set_theme(style="whitegrid")
np.random.seed(SEED)


#---------------------------------HELPERS---------------------------------------

def scatter_plot(df, x, y, title, xlabel=None, ylabel=None):
    """Scatter plot coloured by cell line, shaped by knockdown status."""
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(data=df, x=x, y=y, hue="cell_line", style="knockdown",
                    s=90, alpha=0.85, ax=ax)
    ax.set_title(title)
    ax.set_xlabel(xlabel or x)
    ax.set_ylabel(ylabel or y)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    plt.show()


def prepare_features(df, target):
    """Split a dataframe into X / y. Categorical predictors are coded as integers
    (same behaviour as ranger's default 'ignore' for unordered factors)."""
    X = df.drop(columns=target).copy()
    for col in X.columns[X.dtypes == "category"]:
        X[col] = X[col].cat.codes
    return X, df[target]


def stratified_split(X, y, seed=SEED):
    """75% train / 25% test, stratified on the target (rsample::initial_split)."""
    return train_test_split(X, y, test_size=0.25, stratify=y, random_state=seed)


def oob_permutation_importance(rf, X, y_idx, seed=SEED):
    """Out-of-bag permutation importance (the same idea as ranger's
    importance = 'permutation'): for each tree, the drop in OOB accuracy
    after shuffling a feature, averaged over all trees."""
    rng = np.random.default_rng(seed)
    n, p = X.shape
    all_idx = np.arange(n)
    importance = np.zeros(p)
    inbag_sets = rf.estimators_samples_
    for tree, inbag in zip(rf.estimators_, inbag_sets):
        oob = np.setdiff1d(all_idx, inbag)
        if oob.size == 0:
            continue
        X_oob, y_oob = X[oob], y_idx[oob]
        base_acc = np.mean(tree.predict(X_oob) == y_oob)
        feats = tree.tree_.feature
        # a feature never used in a tree has an importance of exactly 0
        for j in np.unique(feats[feats >= 0]):
            X_perm = X_oob.copy()
            X_perm[:, j] = rng.permutation(X_perm[:, j])
            importance[j] += base_acc - np.mean(tree.predict(X_perm) == y_oob)
    return importance / len(rf.estimators_)


def plot_importance(importance, title, n=20):
    """Horizontal bar plot of the top n variables (like vip::vip)."""
    top = importance.sort_values(ascending=False).head(n)[::-1]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.barh(top.index, top.values, color="grey")
    ax.set_title(title)
    ax.set_xlabel("Importance")
    fig.tight_layout()
    plt.show()


def run_random_forest(rf_df, target, title):
    """Fit a Random Forest (700 trees) on a 75/25 stratified split, then plot
    the confusion matrix, print performance metrics and plot the importance."""
    X, y = prepare_features(rf_df, target)
    classes = list(y.cat.categories)
    X_train, X_test, y_train, y_test = stratified_split(X, y)

    # step_zv(): remove zero-variance predictors (learned on the training set)
    keep = X_train.columns[X_train.var(ddof=0) > 0]
    X_train, X_test = X_train[keep], X_test[keep]

    # mtry = floor(sqrt(ncol(train) - 1)), ncol(train) includes the target
    mtry = int(np.floor(np.sqrt(rf_df.shape[1] - 1)))

    Xtr = X_train.to_numpy(dtype=np.float32)
    Xte = X_test.to_numpy(dtype=np.float32)
    ytr = np.array([classes.index(c) for c in y_train])

    rf = RandomForestClassifier(n_estimators=700, max_features=mtry,
                                random_state=SEED, n_jobs=-1)
    rf.fit(Xtr, ytr)

    # predictions
    pred = np.array(classes, dtype=object)[rf.predict(Xte)]

    # confusion matrix --> plotted as a heatmap
    fig, ax = plt.subplots(figsize=(6, 5))
    ConfusionMatrixDisplay.from_predictions(y_test.astype(str), pred,
                                            labels=classes, cmap="Blues", ax=ax)
    ax.set_title(f"Confusion matrix — {title}")
    fig.tight_layout()
    plt.show()

    # model performance (yardstick::metrics -> accuracy & kappa)
    print(f"\n[{title}] accuracy = {accuracy_score(y_test.astype(str), pred):.3f}, "
          f"kappa = {cohen_kappa_score(y_test.astype(str), pred):.3f}")

    # variable importance (OOB permutation importance)
    importance = pd.Series(oob_permutation_importance(rf, Xtr, ytr),
                           index=keep, name="importance")
    plot_importance(importance, f"Variable importance — {title}", n=20)


def agglomerative_coefficient(Z, n_obs):
    """Agglomerative coefficient (cluster::agnes$ac) from a scipy linkage
    matrix: mean over observations of 1 - (height of first merge / final height)."""
    first_merge = np.zeros(n_obs)
    for a, b, height, _ in Z:
        for node in (int(a), int(b)):
            if node < n_obs:
                first_merge[node] = height
    return float(np.mean(1 - first_merge / Z[-1, 2]))


def hierarchical_clustering(data, labels, label_font_size=6):
    """Compare linkage methods with the agglomerative coefficient, then
    draw the dendrogram obtained with Ward's method."""
    values = data.to_numpy(dtype=float)
    methods = ["average", "single", "complete", "ward"]

    # the closer this value is to 1, the stronger the clusters
    coefs = pd.Series(
        {m: agglomerative_coefficient(linkage(values, method=m), values.shape[0])
         for m in methods})
    print(coefs)

    # hierarchical clustering using the chosen method (here ward is better)
    Z = linkage(values, method="ward")
    fig, ax = plt.subplots(figsize=(14, 6))
    dendrogram(Z, labels=list(labels), leaf_rotation=90,
               leaf_font_size=label_font_size, ax=ax)
    ax.set_title("Dendrogram")
    fig.tight_layout()
    plt.show()


#--------------------------------DATA LOADING-----------------------------------

# import dataset from a .csv file (read.csv2 = ';' separator, ',' decimal)
dataset = "GSE266566_COUNTS.csv"
df_BRCA = pd.read_csv(dataset, sep=";", dtype=str)

#--------------------------------CLEANING---------------------------------------

# rename columns from column 1,2... to actual names (first data row)
df_BRCA.columns = df_BRCA.iloc[0].astype(str).to_list()
df_BRCA = df_BRCA.iloc[1:].reset_index(drop=True)

# delete useless ensemble IDs
df_BRCA = df_BRCA.drop(columns=["Ensembl.114.Transcript.ID"], errors="ignore")

# convert columns from str to numeric to use them (columns 3 to 54 in R)
sample_cols = df_BRCA.columns[2:54]
df_BRCA[sample_cols] = df_BRCA[sample_cols].apply(
    lambda col: pd.to_numeric(col.str.replace(",", ".", regex=False), errors="coerce"))

#-----------------------------------LABELS--------------------------------------

# 52 samples (use the names of the columns as sample names)
sample_names = list(sample_cols)

metadata = pd.DataFrame({
    "sample": sample_names,
    "cell_line": [re.sub(r"counts\.([^.]+)\..*", r"\1", s) for s in sample_names],
    "condition": [re.sub(r"counts\.[^.]+\.([^.]+\.[^.]+)\.RNA.*", r"\1", s)
                  for s in sample_names],
}, index=sample_names)

metadata["knockdown"] = np.where(metadata["condition"].str.contains("NRP1"), "KD", "Control")
metadata["tech"] = np.where(metadata["condition"].str.startswith("sh"), "shRNA", "siRNA")

#-------------------------------NORMALIZATION-----------------------------------

# select only necessary variables (omit name variables), and use the name of
# the transcript as index
expr_mat = df_BRCA[sample_cols].astype(float)
expr_mat.index = df_BRCA["Ensembl.114.Transcript.Name"]
if expr_mat.index.has_duplicates:
    warnings.warn("Duplicated transcript names found: selection by name may be ambiguous.")

# keep transcripts that are present in AT LEAST 3 samples
keep = (expr_mat > 1).sum(axis=1) >= 3
expr_mat = expr_mat[keep]
print("Transcrits retenus après filtre:", expr_mat.shape[0])

# normalization w/ log2(CPM + 1)
lib_sizes = expr_mat.sum(axis=0)
cpm_mat = expr_mat.div(lib_sizes, axis=1) * 1e6
log_cpm = np.log2(cpm_mat + 1)

# automatically adjust n_top_rows parameter
# if number of transcripts post filtration is lower than n_top_rows, use the
# number of transcripts instead
if n_top_rows > log_cpm.shape[0]:
    warnings.warn(f"n_top_rows > nombre de transcrits filtrés. "
                  f"Utilisation de {log_cpm.shape[0]} transcrits à la place.")
    n_top_rows = log_cpm.shape[0]

# selection of the n_top_rows transcripts with the highest variance
variances = log_cpm.var(axis=1, ddof=1).to_numpy()
top_idx = np.argsort(-variances, kind="stable")[:n_top_rows]
log_top = log_cpm.iloc[top_idx]   # FINAL DF TO BE USED (transcripts x samples)

del expr_mat, keep, lib_sizes, cpm_mat, variances, log_cpm, top_idx

#--------------------------------FULL DATA PCA----------------------------------

# samples x transcripts, centered and scaled (prcomp(center = TRUE, scale. = TRUE))
X_samples = log_top.T
X_scaled = (X_samples - X_samples.mean()) / X_samples.std(ddof=1)

pca = PCA(random_state=SEED)
pca_scores = pca.fit_transform(X_scaled.to_numpy())
pc_names = [f"PC{i + 1}" for i in range(pca_scores.shape[1])]

# PCA results on the 8 first dimensions
pca_data = pd.concat(
    [pd.DataFrame(pca_scores[:, :8], columns=pc_names[:8], index=X_samples.index),
     metadata], axis=1)
pca_var_explained = np.round(100 * pca.explained_variance_ratio_, 1)

# PCA plot on PC1 & PC2
scatter_plot(pca_data, "PC1", "PC2", "PCA — NRP1 knockdown vs Control",
             f"PC1 ({pca_var_explained[0]}%)", f"PC2 ({pca_var_explained[1]}%)")

# PCA plot on PC3 & PC4
scatter_plot(pca_data, "PC3", "PC4", "PCA — NRP1 knockdown vs Control",
             f"PC3 ({pca_var_explained[2]}%)", f"PC4 ({pca_var_explained[3]}%)")

# eigenvalues table (eigenvalue = variance of a dimension; also in percent)
pca_eig_val = pd.DataFrame({
    "eigenvalue": pca.explained_variance_,
    "variance.percent": 100 * pca.explained_variance_ratio_,
    "cumulative.variance.percent": 100 * np.cumsum(pca.explained_variance_ratio_),
}, index=[f"Dim.{i + 1}" for i in range(len(pc_names))])

# bar plot of the eigenvalues (scree plot, first 10 dimensions)
fig, ax = plt.subplots(figsize=(8, 5))
scree = pca_eig_val["variance.percent"].head(10)
bars = ax.bar(range(1, len(scree) + 1), scree.values, color="steelblue")
ax.bar_label(bars, labels=[f"{v:.1f}" for v in scree.values], padding=2)
ax.set_ylim(0, 50)
ax.set_xticks(range(1, len(scree) + 1))
ax.set_xlabel("Dimensions")
ax.set_ylabel("Percentage of explained variances")
ax.set_title("Scree plot")
fig.tight_layout()
plt.show()
# at dimension 6, we have 81.6% of the variance explained == keep 6 dimensions for machine learning

# loadings (transcripts x PCs), i.e. R's pca_result$rotation
loadings = pd.DataFrame(pca.components_.T, index=log_top.index, columns=pc_names)

# most contributing descriptors in PC1 (first dim): contribution = 100 * loading^2
contrib_pc1 = (100 * loadings["PC1"] ** 2).sort_values(ascending=False).head(40)
fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(contrib_pc1.index, contrib_pc1.values, color="steelblue")
ax.axhline(100 / loadings.shape[0], color="red", linestyle="--")   # expected average
ax.set_ylabel("Contributions (%)")
ax.set_title("Contribution of variables to Dim-1")
ax.tick_params(axis="x", rotation=90, labelsize=7)
fig.tight_layout()
plt.show()

# sort drivers from PC1 and PC2 (by absolute loading)
pca_loadings = loadings[["PC1", "PC2"]].rename_axis("Variable").reset_index()
pc1_loadings_sorted = pca_loadings.reindex(pca_loadings["PC1"].abs().sort_values(ascending=False).index)
pc2_loadings_sorted = pca_loadings.reindex(pca_loadings["PC2"].abs().sort_values(ascending=False).index)
print("\nTop 10 drivers of PC1:\n", pc1_loadings_sorted.head(10).to_string(index=False))
print("\nTop 10 drivers of PC2:\n", pc2_loadings_sorted.head(10).to_string(index=False))

# remove variables that are not useful for other methods
del pca_data, pca_var_explained, pca_eig_val, pca_loadings, pc1_loadings_sorted, pc2_loadings_sorted

#------------------------------------UMAP---------------------------------------

# UMAP on the (unscaled) samples x transcripts matrix
umap_layout = umap.UMAP(n_components=2, random_state=SEED).fit_transform(
    X_samples.to_numpy())
umap_df = pd.concat(
    [pd.DataFrame(umap_layout, columns=["UMAP1", "UMAP2"], index=X_samples.index),
     metadata], axis=1)

# UMAP plot
scatter_plot(umap_df, "UMAP1", "UMAP2", "UMAP — NRP1 knockdown vs Control")

del umap_layout, umap_df

#-----------------------------PCA LOADINGS EXTRACT------------------------------

# choose Principal Components to consider, we take the first 6 --> see eigenvalues
pcs_to_use = pc_names[:6]

# compute importance score per gene
gene_scores = loadings[pcs_to_use].abs().sum(axis=1)

# select top genes
top_gene_n = min(500, len(gene_scores))
top_genes = gene_scores.sort_values(ascending=False, kind="stable").index[:top_gene_n]

# dataframe with only the top genes from log_top, genes as columns
rf_df = log_top.loc[top_genes].T.copy()

del loadings, pcs_to_use, gene_scores, top_gene_n

#------------------ML: CLASSIFICATION USING RANDOM FOREST-----------------------

# add knockdown target col and cell line col from metadata
rf_df["knockdown"] = pd.Categorical(metadata.loc[rf_df.index, "knockdown"],
                                    categories=["Control", "KD"])
rf_df["cell_line"] = pd.Categorical(metadata.loc[rf_df.index, "cell_line"])

# Random Forest with permutation importance: measures how much the model's
# prediction error increases when the values of a feature are shuffled
# == helps identify which features are most important for the predictions

#-----------------Knockdown recognition model
run_random_forest(rf_df, target="knockdown", title="Knockdown recognition")

#------Cell_line recognition model
run_random_forest(rf_df, target="cell_line", title="Cell line recognition")

#--------------------------HIERARCHICAL CLUSTERING------------------------------

#--------cell lines using genes (samples are clustered)

# drop cell_line and knockdown columns --> numeric matrix only
df_clust_cl = rf_df.drop(columns=["knockdown", "cell_line"])

# scale dataframe since the genes have different variance ranges
df_clust_cl = (df_clust_cl - df_clust_cl.mean()) / df_clust_cl.std(ddof=1)

hierarchical_clustering(df_clust_cl, labels=df_clust_cl.index, label_font_size=8)

del df_clust_cl

#--------genes using cell lines (genes are clustered)

# scale each sample column across the genes
df_clust_kd = log_top.loc[top_genes]
df_clust_kd = (df_clust_kd - df_clust_kd.mean()) / df_clust_kd.std(ddof=1)

hierarchical_clustering(df_clust_kd, labels=df_clust_kd.index, label_font_size=4)

del df_clust_kd

#------------------------------------C5.0---------------------------------------
# NOTE: there is no maintained C5.0 implementation in scikit-learn. The closest
# equivalent is a CART decision tree using the entropy (information gain)
# criterion, which is what is used here.

X_c50, y_c50 = prepare_features(rf_df, target="cell_line")
X_train, X_test, y_train, y_test = stratified_split(X_c50, y_c50)

model = DecisionTreeClassifier(criterion="entropy", random_state=SEED)
model.fit(X_train, y_train.astype(str))

# textual summary of the tree (to check the column names)
print("\n" + export_text(model, feature_names=list(X_train.columns)))
print(f"Test accuracy: {model.score(X_test, y_test.astype(str)):.3f}")

# decision tree graph
fig, ax = plt.subplots(figsize=(16, 8))
plot_tree(model, feature_names=list(X_train.columns), class_names=list(model.classes_),
          filled=True, rounded=True, fontsize=7, ax=ax)
fig.tight_layout()
plt.show()
