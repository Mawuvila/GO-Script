### Gene enrichment analysis 
### Olink data
### Hydrops 
### Date = 10 May 2024


### Install pakages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("treeio")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationDbi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

library(treeio)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)


setwd("path")
# load the deseq_results.csv data
Data <- read.csv("data.csv", header = TRUE, sep = ",", row.names = 1)
head(Data)[1-5, ]

## take only DE genes 
# filter for genes with padj < 0.05
sigs_filtered <- sigs %>% filter(sigs$adj_pval < 0.05) 


# filter for genes with log2fold < 1 and less than -1
#sigs_filtered <- sigs_filtered %>% filter(abs(sigs_filtered$log2FoldChange) > 0.5)

## get the row names 
genes_to_test <- rownames(sigs_filtered)
genes_to_test

### Run the go command. >>  ont should be One of "MF", "BP", and "CC" subontologies.
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
head(GO_results)


## convert the GO_results into a data frame 
as.data.frame(GO_results)
write.csv(GO_results, "GO_results.csv")
head(GO_results)[1, ]


### bar plot
fit <- plot(barplot(GO_results, showCategory = 20))
png("out.png", res = 250, width = 1200, height = 1000)
fit


### Volcano plots
## install package
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)


### bubble plot
## install necessary package 
install.packages("tidyverse")
install.packages("plotly")
library(tidyverse)
library(plotly)
# plot the bubble plot
Bplot <- GO_results %>% ggplot(aes(x = Count, y = Description, size = p.adjust))
Bplot





### run volcano plot with data loaded here 
deseq_results <- read.csv("deseq_results.csv", header = TRUE, sep = ",")
EnhancedVolcano(deseq_results, x = "log2FoldChange", y = "padj", lab = deseq_results$X)

#change the thresholds 
EnhancedVolcano(deseq_results, x = "log2FoldChange", y = "padj", lab = deseq_results$X, 
                pCutoff = 1e-1, FCcutoff = 0.5)

#lable only some selected genes
selected_genes <- c("IGLV3-10", "IGHV6-1", "IGKV1D-13", "APOE", "IGKV1-16", "IGKV1D-39", "PZP", "LRP1")
EnhancedVolcano(deseq_results, x = "log2FoldChange", y = "padj", lab = deseq_results$X, 
                pCutoff = 1e-1, FCcutoff = 1, 
                selectLab = selected_genes)



