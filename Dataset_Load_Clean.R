
###Loading the dataset###

df_labels <- read.csv2("labels.csv", sep = ",")
row_list <- df_labels$X
df_data <- read.csv2("data.csv", sep = ",",row.names = row_list)
df_data$X <- NULL
df_final <- cbind(df_labels$Class,df_data)  
names(df_final)[names(df_final) == "df_labels$Class"] <-"Tumor_type"

df_BRCA <- df_final[df_final$Tumor_type == 'BRCA',]

write.csv(df_BRCA,"BRCA_df.csv",
          quote=FALSE,sep = " ",eol = "\r\n",col.names = names(df_BRCA)
          )

df_test <-read.csv2("BRCA_only.csv", sep = ",")
df_test2 <-read.csv2("BRCA_only_2.csv", sep = ",",row.names = )
write.