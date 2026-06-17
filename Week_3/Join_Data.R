library(dplyr)

#---------------------------------PARAMETERS------------------------------------

n_top_rows <- 2000   #numbers of transcripts to be used

df_BRCA <- read.csv2("../Datasets/GSE266566_full.csv")
df_HD <-read.csv2("GSE71862_MCF7_MCF10A_RSEM_expectedcounts.csv")

#--------------------------------CLEANING---------------------------------------

#rename columns from column 1,2... to actual names
colnames(df_BRCA) <- df_BRCA[1,] 
df_BRCA <- df_BRCA[-1,]
df_BRCA <- df_BRCA[ ,c(-1,-3)]


colnames(df_HD) <- df_HD[1,]
df_HD <- df_HD[-1,]

#delete useless columns
df_HD$accession <- NULL

#convert column from chr to numeric to use them
df_BRCA[,2:53] <- lapply(df_BRCA[,2:53], function(x) as.numeric(as.character(x)))

#convert column from chr to numeric to use them
df_HD[,2:7] <- lapply(df_HD[,2:7], function(x) as.numeric(as.character(x)))

#--------------------------TREATMENT AND COMBINATION-----------------------------

df_BRCA$gene <- sub('([^-]+).*','\\1',df_BRCA$Ensembl.114.Transcript.Name)
df_BRCA$total <-  rowSums(df_BRCA[,c(2:53)])

df_BRCA_reduced <- slice_max(df_BRCA, order_by= total , by=gene)
df_BRCA_reduced$Ensembl.114.Transcript.Name <- df_BRCA_reduced$gene
df_BRCA_reduced$gene <-NULL

names(df_BRCA_reduced)[names(df_BRCA_reduced) == 'Ensembl.114.Transcript.Name'] <- 'gene'

df_combined <- inner_join(df_BRCA_reduced,df_HD)
df_combined$total <- NULL

#-----------------------------------LABELS--------------------------------------                   

# 52 samples (use the names of the column as sample names)
sample_names <- colnames(df_combined)[2:59] 

metadata <- data.frame(
  sample = sample_names,
  cell_line = gsub('counts\\.|\\.s.*|_.*',"",sample_names),
  stringsAsFactors = FALSE
)

metadata$knockdown <- ifelse(grepl("MCF10A", metadata$cell_line), "Healthy", "Cancerous")

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


