#---------------------------------LOADING DATA----------------------------------

#Loading libraries
library(dplyr)

#Check loading of libraries
sessionInfo()

#Dataset we will be using 
dataset <- "GSE266566_COUNTS.rsem_human.transcripts.csv"  #path of the dataset

#Loading a dataset from a csv file into a dataframe called df_BRCA
df_BRCA <- read.csv2(dataset)

#Remove variables that are no longer needed
rm(dataset)

#Rename columns 
head(df_BRCA)
colnames(df_BRCA) <- df_BRCA[c(1),]
df_BRCA <- df_BRCA[-c(1),]

#Delete Ensembl ID
df_BRCA$Ensembl.114.Transcript.ID <- NULL

#Check for missing values
values <- colSums(is.na(df_BRCA))
check_missing <- data.frame(values)

sample_names <- row.names(check_missing)
for(i in range(length(check_missing$values))){
  print(check_missing$values[i])
  if(check_missing$values[i] !=0){
    print("Value missing in:",sample_names[i])
  }else {print("There are no missing values")
  }
}

#Remove variables that are no longer needed
rm(check_missing, i, sample_names, values)

#Convert chr cols to num type for PCA
df_BRCA[, 3:54] <- as.numeric(unlist(df_BRCA[, 3:54])) 

#Checking the type of the data in the dataframe
str(df_BRCA)

#Checking the dataset was correctly imported and cleaned
View(df_BRCA)

#-------------------------------------PCA---------------------------------------

#Loading libraries
library(ggplot2)

#Select only necessary variables (omit name variables)
X <- df_BRCA[, 3:54]

#Standardize data
scaled_X <- scale(X)

#Compute PCA
pca_result <- prcomp(scaled_X, center = TRUE, scale. = TRUE)
pca_result

#Extract principal components for each observation
pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data

#Add Description/Name descriptor for classification
pca_data$Name <- df_BRCA$Ensembl.114.Transcript.Name
pca_data

#Create PCA graph using ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = Name)) +
  geom_point(size = 3) +
  labs(title = "Analyse en Composantes Principales (PCA) sur df_BRCA",
       x = "Première Composante Principale (PC1)",
       y = "Deuxième Composante Principale (PC2)") +
  theme_minimal()

#Remove variables that are no longer needed
rm(pca_data,pca_result,scaled_X,X)




