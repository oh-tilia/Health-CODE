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


#-----------------------------KEEP TOP VARIABLE GENES---------------------------


