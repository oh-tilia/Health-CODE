#-----Package installation-----
install.packages("BiocManager")
install.packages("RColorBrewer")
install.packages("pheatmap")
install.packages("ggrepel")
install.packages("devtools")
install.packages("tidyverse")
install.packages("pandas")

BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("DOSE")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")
BiocManager::install("DEGreport")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("AnnotationHub")
BiocManager::install("ensembldb")

library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(DEGreport)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(ensembldb)
library(pandas)
sessionInfo()

#-----Actual Code-----


countData <- read.csv('GSE266566_Test.csv', header = TRUE, sep = ";")
countData
head(countData)
colnames(countData) <- countData[c(1),]
countData <- countData[-c(1),]
rownames(countData) <- countData$Ensembl.114.Transcript.Name
countData$Ensembl.114.Transcript.Name <- NULL

metaData <- read.csv("GSE266566_Metadata.csv",sep=";")

dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~Ensembl.114.Transcript.Type)
