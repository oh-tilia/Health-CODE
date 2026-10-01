# Healthy (epithelial) cells Data base : https://cellxgene.cziscience.com/collections/c9706a92-0e5f-46c1-96d8-20e42467f287
# Cancer cells Data base : https://cellxgene.cziscience.com/collections/dea97145-f712-431c-a223-6b5f565f362a

library(SingleCellExperiment)
library(SummarizedExperiment)
library(Matrix)
library(anndataR)
df <- read_h5ad(
  path = "scRNAseq_cancer.h5ad",
  as = "SingleCellExperiment",
  mode = "r"
)

genes <- data.frame(rownames(df),rowData(df)$feature_name)


write.table(genes,
            file = "genes.tsv",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(colnames(df),
            file = "barcodes.tsv",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# extraire la matrice
mat <- assay(df,"X")
mat <- as(mat, "dgCMatrix")
writeMM(mat, file = "matrix.mtx")

# Rendre les noms de gènes uniques AVANT CreateSeuratObject
# Récupérer les symboles
genes <- rowData(df)$feature_name

# Convertir en caractères
genes <- as.character(genes)

# Remplacer les NA / vides par les IDs Ensembl
genes[is.na(genes) | genes == ""] <- rownames(df)[is.na(genes) | genes == ""]

# Rendre les noms uniques
genes <- make.unique(genes)

# Appliquer les noms à la matrice
rownames(mat) <- genes

#Créer l'objet Seurat propre
library(Seurat)
seurat_obj <- CreateSeuratObject(
  counts = mat,
  meta.data = as.data.frame(colData(df))
)

seurat_obj
head(seurat_obj@meta.data)

# Normalisation
seurat_obj <- NormalizeData(seurat_obj)

# HVGs (gènes les plus variables)
seurat_obj <- FindVariableFeatures(seurat_obj)

# Scaling
seurat_obj <- ScaleData(seurat_obj)

# Vérifier que la PCA existe
seurat_obj@reductions$pca

# PCA
seurat_obj <- RunPCA(seurat_obj)

# UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

DimPlot(seurat_obj, reduction = "umap", group.by = "celltype_major", label = TRUE)

