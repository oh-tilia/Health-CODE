#==============================================================================#
#                        PIPELINE Cancerous vs Healthy                         #
#==============================================================================#
# Python 3.12 translation of Join_Data.R
#
# Requirements:
#   pip install numpy pandas scikit-learn matplotlib seaborn
#
# All figures are collected while the pipeline runs and shown at the end in ONE
# window: use the <- / -> arrow keys to move between figures (Home / End jump
# to the first / last one). No need to close a figure to see the next one.

#---------------------------------LIBRARIES-------------------------------------

import re
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier, export_text, plot_tree
from sklearn.model_selection import train_test_split
from sklearn.metrics import (accuracy_score, cohen_kappa_score,
                             ConfusionMatrixDisplay)

sns.set_theme(style="whitegrid")

#---------------------------------PARAMETERS------------------------------------

n_top_rows = 2000            # numbers of transcripts to be used
SEED = 222                   # seed for the Random Forest and the first C5.0 model
SEED_C50_STATE = 666         # seed for the last C5.0 model (state recognition)
SHOW_PC1_CONTRIBUTIONS = False   # optional PCA drivers (commented out in the R code)

BRCA_FILE = "GSE266566_COUNTS.csv"
HD_FILE = "GSE71862_MCF7_MCF10A_RSEM_expectedcounts.csv"


#-------------------------FIGURE BROWSER (arrow navigation)---------------------

class FigureBrowser:
    """Collects plots (as drawing functions) and displays them in a single
    window. <- / -> switch figure, Home / End go to the first / last one."""

    def __init__(self):
        self.plots = []      # list of (name, draw_function)
        self.index = 0
        self.fig = None

    def add(self, draw, name):
        """Register a plot: `draw(fig)` must draw on the figure it receives."""
        self.plots.append((name, draw))

    def render(self):
        name, draw = self.plots[self.index]
        fig = self.fig
        fig.clf()
        draw(fig)
        fig.text(0.5, 0.008,
                 f"{self.index + 1}/{len(self.plots)} - {name}     "
                 "(<- / -> : previous / next,  Home / End : first / last)",
                 ha="center", va="bottom", fontsize=9, color="dimgray")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig.tight_layout(rect=(0, 0.035, 1, 1))
        manager = fig.canvas.manager
        if manager is not None:
            manager.set_window_title(f"Figure {self.index + 1}/{len(self.plots)} - {name}")
        fig.canvas.draw_idle()

    def on_key(self, event):
        last = len(self.plots) - 1
        if event.key == "right":
            self.index = min(self.index + 1, last)
        elif event.key == "left":
            self.index = max(self.index - 1, 0)
        elif event.key == "home":
            self.index = 0
        elif event.key == "end":
            self.index = last
        else:
            return
        self.render()

    def show(self):
        if not self.plots:
            return
        # free the arrow / home keys from matplotlib's default toolbar shortcuts
        for param, key in (("keymap.back", "left"), ("keymap.forward", "right"),
                           ("keymap.home", "home")):
            plt.rcParams[param] = [k for k in plt.rcParams[param] if k != key]
        self.fig = plt.figure(figsize=(12, 7.5))
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.render()
        plt.show()


browser = FigureBrowser()


#---------------------------------PLOT HELPERS----------------------------------

def scatter_plot(df, x, y, title, xlabel=None, ylabel=None, name=None):
    """Scatter plot coloured by cell line, shaped by state (cancerous/healthy)."""
    def draw(fig):
        ax = fig.add_subplot()
        sns.scatterplot(data=df, x=x, y=y, hue="cell_line", style="state",
                        s=90, alpha=0.85, ax=ax)
        ax.set_title(title)
        ax.set_xlabel(xlabel or x)
        ax.set_ylabel(ylabel or y)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    browser.add(draw, name or title)


def plot_scree(var_percent, n_dims=10):
    """Bar plot of the percentage of variance explained by each dimension."""
    values = np.asarray(var_percent)[:n_dims]

    def draw(fig):
        ax = fig.add_subplot()
        bars = ax.bar(range(1, len(values) + 1), values, color="steelblue")
        ax.bar_label(bars, labels=[f"{v:.1f}" for v in values], padding=2)
        ax.set_ylim(0, 50)
        ax.set_xticks(range(1, len(values) + 1))
        ax.set_xlabel("Dimensions")
        ax.set_ylabel("Percentage of explained variances")
        ax.set_title("Scree plot")
    browser.add(draw, "Scree plot")


def plot_contributions(loadings, pc="PC1", top=40):
    """Contribution (%) of the top variables to one principal component."""
    contrib = (100 * loadings[pc] ** 2).sort_values(ascending=False).head(top)
    expected = 100 / loadings.shape[0]   # contribution if all variables were equal

    def draw(fig):
        ax = fig.add_subplot()
        ax.bar(range(len(contrib)), contrib.values, color="steelblue")
        ax.axhline(expected, color="red", linestyle="--")
        ax.set_xticks(range(len(contrib)))
        ax.set_xticklabels(contrib.index, rotation=90, fontsize=7)
        ax.set_ylabel("Contributions (%)")
        ax.set_title(f"Contribution of variables to {pc}")
    browser.add(draw, f"{pc} contributions")


def plot_confusion(y_true, y_pred, labels, title):
    """Confusion matrix drawn as a heatmap."""
    def draw(fig):
        ax = fig.add_subplot()
        ConfusionMatrixDisplay.from_predictions(y_true, y_pred, labels=labels,
                                                cmap="Blues", ax=ax)
        ax.grid(False)
        ax.set_title(f"Confusion matrix - {title}")
    browser.add(draw, f"Confusion matrix - {title}")


def plot_importance(importance, title, n=20):
    """Horizontal bar plot of the top n variables (like vip::vip)."""
    top = importance.sort_values(ascending=False).head(n)[::-1]

    def draw(fig):
        ax = fig.add_subplot()
        ax.barh(list(top.index), top.values, color="grey")
        ax.set_xlabel("Importance")
        ax.set_title(f"Variable importance - {title}")
    browser.add(draw, f"Variable importance - {title}")


def plot_decision_tree(model, feature_names, title):
    """Decision tree graph."""
    def draw(fig):
        ax = fig.add_subplot()
        plot_tree(model, feature_names=feature_names, class_names=list(model.classes_),
                  filled=True, rounded=True, fontsize=7, ax=ax)
        ax.set_title(title)
    browser.add(draw, title)


#---------------------------------DATA HELPERS----------------------------------

def r_make_names(names):
    """Approximation of R's make.names(unique = TRUE), which read.csv2 applies
    to the header of a file."""
    fixed = []
    for name in names:
        name = "" if pd.isna(name) else str(name)
        name = re.sub(r"[^0-9A-Za-z._]", ".", name)
        if not re.match(r"([A-Za-z]|\.(?![0-9]))", name):
            name = "X" + name
        fixed.append(name)
    seen, unique = {}, []
    for name in fixed:
        if name in seen:
            seen[name] += 1
            unique.append(f"{name}.{seen[name]}")
        else:
            seen[name] = 0
            unique.append(name)
    return unique


def read_csv2(path):
    """read.csv2 (';' separator, ',' decimal) + the 'Column1' header fix:
    when the header is generic (Column1, Column2...), the first data row holds
    the real column names."""
    raw = pd.read_csv(path, sep=";", dtype=str, header=None)
    columns, data = r_make_names(raw.iloc[0]), raw.iloc[1:]
    if columns[0] == "Column1":
        columns, data = raw.iloc[1].astype(str).to_list(), raw.iloc[2:]
    data = data.copy()
    data.columns = columns
    return data.reset_index(drop=True)


def numeric_columns(df, start, stop):
    """as.numeric(as.character(x)) on df.iloc[:, start:stop] (decimal commas accepted)."""
    converted = df.iloc[:, start:stop].apply(
        lambda col: pd.to_numeric(col.str.replace(",", ".", regex=False), errors="coerce"))
    return pd.concat([df.iloc[:, :start], converted, df.iloc[:, stop:]], axis=1)


#---------------------------------ML HELPERS------------------------------------

def prepare_features(df, target):
    """Split a dataframe into X / y. Categorical predictors are coded as integers
    (same behaviour as ranger's default 'ignore' for unordered factors)."""
    X = df.drop(columns=target).copy()
    for col in X.columns[X.dtypes == "category"]:
        X[col] = X[col].cat.codes
    return X, df[target]


def stratified_split(X, y, seed):
    """75% train / 25% test, stratified on the target (rsample::initial_split)."""
    return train_test_split(X, y, test_size=0.25, stratify=y, random_state=seed)


def oob_permutation_importance(rf, X, y_idx, seed):
    """Out-of-bag permutation importance (the same idea as ranger's
    importance = 'permutation'): for each tree, the drop in OOB accuracy
    after shuffling a feature, averaged over all trees."""
    rng = np.random.default_rng(seed)
    n, p = X.shape
    all_idx = np.arange(n)
    importance = np.zeros(p)
    for tree, inbag in zip(rf.estimators_, rf.estimators_samples_):
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


def run_random_forest(rf_df, target, title, seed):
    """Fit a Random Forest (700 trees) on a 75/25 stratified split, print a
    summary, then plot the confusion matrix and the variable importance."""
    X, y = prepare_features(rf_df, target)
    classes = list(y.cat.categories)
    X_train, X_test, y_train, y_test = stratified_split(X, y, seed)

    # step_zv(): remove zero-variance predictors (learned on the training set)
    keep = X_train.columns[X_train.var(ddof=0) > 0]
    X_train, X_test = X_train[keep], X_test[keep]

    # mtry = floor(sqrt(ncol(train) - 1)), ncol(train) includes the target
    mtry = int(np.floor(np.sqrt(rf_df.shape[1] - 1)))

    Xtr = X_train.to_numpy(dtype=np.float32)
    Xte = X_test.to_numpy(dtype=np.float32)
    ytr = y_train.cat.codes.to_numpy()

    rf = RandomForestClassifier(n_estimators=700, max_features=mtry, oob_score=True,
                                random_state=seed, n_jobs=-1)
    rf.fit(Xtr, ytr)

    # summary of the underlying model (like printing the ranger object)
    print(f"\n[{title}] Random Forest")
    print(f"  Type:                       Classification")
    print(f"  Number of trees:            {rf.n_estimators}")
    print(f"  Sample size:                {Xtr.shape[0]}")
    print(f"  Number of independent vars: {Xtr.shape[1]}")
    print(f"  Mtry:                       {mtry}")
    print(f"  Target node size:           {rf.min_samples_leaf}")
    print(f"  Variable importance mode:   permutation")
    print(f"  Splitrule:                  {rf.criterion}")
    print(f"  OOB prediction error:       {100 * (1 - rf.oob_score_):.2f} %")

    # predictions
    pred = np.array(classes, dtype=object)[rf.predict(Xte)]
    y_true = y_test.astype(str)

    # confusion matrix --> plotted as a heatmap
    plot_confusion(y_true, pred, classes, title)

    # model performance (yardstick::metrics -> accuracy & kappa)
    print(f"[{title}] accuracy = {accuracy_score(y_true, pred):.3f}, "
          f"kappa = {cohen_kappa_score(y_true, pred):.3f}")

    # variable importance (OOB permutation importance)
    importance = pd.Series(oob_permutation_importance(rf, Xtr, ytr, seed),
                           index=keep, name="importance")
    plot_importance(importance, title, n=20)


def run_decision_tree(rf_df, target, title, seed):
    """Decision tree on a 75/25 stratified split. There is no maintained C5.0
    implementation in scikit-learn: the closest equivalent is a CART tree using
    the entropy (information gain) criterion, which is what is used here."""
    X, y = prepare_features(rf_df, target)
    X_train, X_test, y_train, y_test = stratified_split(X, y, seed)

    model = DecisionTreeClassifier(criterion="entropy", random_state=seed)
    model.fit(X_train, y_train.astype(str))

    # text summary of the tree --> verify genes used
    print(f"\n[{title}] Decision tree")
    print(export_text(model, feature_names=list(X_train.columns), max_depth=50))
    print(f"Training accuracy: {model.score(X_train, y_train.astype(str)):.3f}")
    print(f"Test accuracy:     {model.score(X_test, y_test.astype(str)):.3f}")

    # decision tree graph
    plot_decision_tree(model, list(X_train.columns), title)


#--------------------------------DATA LOADING-----------------------------------

df_BRCA = read_csv2(BRCA_FILE)
df_HD = read_csv2(HD_FILE)

#--------------------------------CLEANING---------------------------------------

# remove unneeded columns from the dataset (columns 1 and 3 of BRCA, 'accession' of HD)
df_BRCA = df_BRCA.iloc[:, np.delete(np.arange(df_BRCA.shape[1]), [0, 2])]
df_HD = df_HD.drop(columns="accession", errors="ignore")

# convert columns from str to numeric in order to use them
df_BRCA = numeric_columns(df_BRCA, 1, 53)
df_HD = numeric_columns(df_HD, 1, 7)

#--------------------------TREATMENT AND COMBINATION-----------------------------

# regex to get only the gene name from the transcript name (drop the "-xxx" suffix)
df_BRCA["gene"] = df_BRCA["Ensembl.114.Transcript.Name"].str.replace(
    r"([^-]+).*", r"\1", n=1, regex=True)

# total expression of each transcript across all samples
# --> if duplicate genes, only the one with the most expression is kept (slice_max)
df_BRCA["total"] = df_BRCA.iloc[:, 1:53].sum(axis=1)

is_max = df_BRCA["total"] == df_BRCA.groupby("gene")["total"].transform("max")
group_order = df_BRCA.groupby("gene", sort=False).ngroup()   # order of first appearance
df_BRCA_reduced = (df_BRCA[is_max]
                   .assign(_grp=group_order[is_max])
                   .sort_values("_grp", kind="stable")
                   .drop(columns="_grp")
                   .copy())

# replace the column with the transcript name by the gene name (kept on the left)
df_BRCA_reduced["Ensembl.114.Transcript.Name"] = df_BRCA_reduced["gene"]
df_BRCA_reduced = (df_BRCA_reduced.drop(columns="gene")
                   .rename(columns={"Ensembl.114.Transcript.Name": "gene"}))

# combine the two datasets: 3 healthy samples and 50+ cancerous samples
common_cols = [c for c in df_BRCA_reduced.columns if c in df_HD.columns]
print("Joining by:", common_cols)
df_combined = df_BRCA_reduced.merge(df_HD, how="inner", on=common_cols).drop(columns="total")

del df_BRCA, df_BRCA_reduced, df_HD, is_max, group_order

#-----------------------------------LABELS--------------------------------------

# 58 samples (use the names of the columns as sample names)
sample_names = list(df_combined.columns[1:59])

# metadata dataframe to easily access the cell line and state of each sample
metadata = pd.DataFrame({
    "sample": sample_names,
    "cell_line": [re.sub(r"counts\.|\.s.*|_.*", "", s) for s in sample_names],
})

# add state (cancerous or healthy) of cells to the metadata
metadata["state"] = np.where(metadata["cell_line"].str.contains("MCF10A"),
                             "Healthy", "Cancerous")

#-------------------------------NORMALIZATION-----------------------------------

# select only necessary variables (omit name variables), genes as row names
expr_mat = df_combined.iloc[:, 1:59].astype(float)
expr_mat.index = df_combined["gene"].to_numpy()
if expr_mat.index.has_duplicates:
    warnings.warn("Duplicated gene names found after the join (kept, as in the R code).")

# keep transcripts that are present in AT LEAST 3 samples
keep = (expr_mat > 1).sum(axis=1) >= 3
expr_mat = expr_mat[keep.to_numpy()]
print("Transcrits retenus après filtre:", expr_mat.shape[0])

# normalization w/ log2(CPM + 1)
# Log2 is used because it aids in calculating fold change (up- vs down-regulated
# genes between samples) and is closer to the biologically-detectable changes.
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
log_top = log_cpm.iloc[top_idx]   # FINAL DF TO BE USED (genes x samples)

del expr_mat, keep, lib_sizes, cpm_mat, variances, log_cpm, top_idx, df_combined

#--------------------------------FULL DATA PCA----------------------------------

# samples x genes, centered and scaled (prcomp(center = TRUE, scale. = TRUE))
X_samples = log_top.to_numpy().T
X_scaled = (X_samples - X_samples.mean(axis=0)) / X_samples.std(axis=0, ddof=1)

pca = PCA(random_state=SEED)
pca_scores = pca.fit_transform(X_scaled)
pc_names = [f"PC{i + 1}" for i in range(pca_scores.shape[1])]

# PCA results on the 8 first dimensions
pca_data = pd.concat(
    [pd.DataFrame(pca_scores[:, :8], columns=pc_names[:8]), metadata], axis=1)
pca_var_explained = np.round(100 * pca.explained_variance_ratio_, 1)

# PCA plot on PC1 & PC2
scatter_plot(pca_data, "PC1", "PC2", "PCA - cancerous vs healthy",
             f"PC1 ({pca_var_explained[0]}%)", f"PC2 ({pca_var_explained[1]}%)",
             name="PCA - PC1 vs PC2")

# PCA plot on PC3 & PC4
scatter_plot(pca_data, "PC3", "PC4", "PCA - cancerous vs healthy",
             f"PC3 ({pca_var_explained[2]}%)", f"PC4 ({pca_var_explained[3]}%)",
             name="PCA - PC3 vs PC4")

# eigenvalues table (eigenvalue = variance of a dimension; also in percent)
pca_eig_val = pd.DataFrame({
    "eigenvalue": pca.explained_variance_,
    "variance.percent": 100 * pca.explained_variance_ratio_,
    "cumulative.variance.percent": 100 * np.cumsum(pca.explained_variance_ratio_),
}, index=[f"Dim.{i + 1}" for i in range(len(pc_names))])

# bar plot of the eigenvalues
plot_scree(pca_eig_val["variance.percent"])
# at dimension 6, we have 88.4% of the variance explained == keep 6 dimensions for machine learning

# loadings (genes x PCs), i.e. R's pca_result$rotation
loadings = pd.DataFrame(pca.components_.T, index=log_top.index, columns=pc_names)

# optional (commented out in the R code): most contributing genes in PC1 and
# the drivers of PC1 / PC2 sorted by absolute loading
if SHOW_PC1_CONTRIBUTIONS:
    plot_contributions(loadings, "PC1", top=40)
    pca_loadings = loadings[["PC1", "PC2"]].rename_axis("Variable").reset_index()
    for pc in ("PC1", "PC2"):
        order = pca_loadings[pc].abs().sort_values(ascending=False).index
        print(f"\nTop 10 drivers of {pc}:\n", pca_loadings.loc[order].head(10).to_string(index=False))

del pca_data, pca_var_explained, pca_eig_val

#-----------------------------PCA LOADINGS EXTRACT------------------------------

# choose Principal Components to consider, we take the first 6 --> see eigenvalues
pcs_to_use = pc_names[:6]

# compute importance score per gene
gene_scores = loadings[pcs_to_use].abs().sum(axis=1)

# select top genes
top_gene_n = min(500, len(gene_scores))
top_genes = gene_scores.sort_values(ascending=False, kind="stable").index[:top_gene_n]

# dataframe with only the top genes from log_top, genes as columns.
# Rows are picked by gene name (first match), as R does with log_top[top_genes, ]
first_row = pd.Series(np.arange(len(log_top)), index=log_top.index)
first_row = first_row[~first_row.index.duplicated()]
rows = log_top.to_numpy()[first_row.loc[top_genes].to_numpy()]     # genes x samples

# removing columns with duplicate data (a gene name that was selected twice)
is_dup = pd.DataFrame(rows).duplicated().to_numpy()
print(f"{is_dup.sum()} duplicated column(s) removed")
rf_df = pd.DataFrame(rows[~is_dup].T, index=log_top.columns, columns=top_genes[~is_dup])

del loadings, pcs_to_use, gene_scores, top_gene_n, first_row, rows, is_dup

#------------------ML: CLASSIFICATION USING RANDOM FOREST-----------------------

# add state target col and cell line col from metadata
rf_df["state"] = pd.Categorical(metadata["state"].to_numpy())
rf_df["cell_line"] = pd.Categorical(metadata["cell_line"].to_numpy())

#-----------------State recognition model (cancerous vs healthy)
run_random_forest(rf_df, target="state", title="State recognition", seed=SEED)

#------------------------------------C5.0---------------------------------------
# (decision tree with the entropy criterion, see run_decision_tree)

#-----------------Cell line recognition model
run_decision_tree(rf_df, target="cell_line", title="Decision tree - cell line", seed=SEED)

#-----------------State recognition model, genes only
# remove cell_line --> if not removed the model uses the cell_line label to
# separate states and does not rely on the genes
rf_df_state = rf_df.drop(columns="cell_line")
run_decision_tree(rf_df_state, target="state", title="Decision tree - state",
                  seed=SEED_C50_STATE)

#------------------------------------FIGURES------------------------------------
# open the single window: <- / -> to navigate between all the figures
browser.show()
