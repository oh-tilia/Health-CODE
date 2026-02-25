#File for loading the dataset

#Get necessary files and set them to variables
label_file = ""
dataset_file = ""

#Load Labels for the dataset
df_labels <- read.csv2(label_file, sep = ",")

#Get the labels
row_list <- df_labels$X

#Load dataset with row names set as labels
df_data <- read.csv2("data.csv", sep = ",",row.names = row_list)
df_data$X <- NULL

#Add the types of tumor to the dataset to sort later
df_final <- cbind(df_labels$Class,df_data)  
names(df_final)[names(df_final) == "df_labels$Class"] <-"Tumor_type"


#Create a dataset with only breast cancer related data
df_BRCA <- df_final[df_final$Tumor_type == 'BRCA',]

#Export dataframe as csv 
output_name = "test"
write.csv(df_BRCA,output_name)

#Import dataset to check if it corectly exported
df_test <- read.csv2(output_name, sep = ",")