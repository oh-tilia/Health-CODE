#Loading libraries
library(dplyr)

#Dataset we will be using 
dataset <- "./Datasets/BRCA_Gene_Reads_Healthy_Cropped.csv"  #path of the dataset
 
#Loading a dataset from a txt file 
#df_BRCA <- read.table(dataset)

#Loading a dataset from a csv file into a dataframe called df_BRCA
df_BRCA <- read.csv2(dataset)



#check for missing value
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

#checks for duplicate data in the name of the gene
anyDuplicated(df_BRCA$Name) 

#checks for duplicate data in the description of the gene
anyDuplicated(df_BRCA$Description)

#Viewing the column names to check the attributes
col_names <- colnames(df_BRCA)

#Checking the dataset was correctly imported
View(df_BRCA)

#transform values from the name and description column from character to factor
#this could help later on when visualizing data using certain methods such as PCA
#columns "Name" and "Description" are specific to our test dataset but need to be adapted
df_BRCA <- df_BRCA %>% mutate_at(c('Name','Description'), as.factor)

#Check correct mutation of data types (change the name after the $ to the actual column names)
class(df_BRCA$Name)
class(df_BRCA$Description)

#Loop to check the classes of all the column in the dataframe
#for(names in colnames(df_BRCA)){print(class(df_BRCA[[names]]))}












