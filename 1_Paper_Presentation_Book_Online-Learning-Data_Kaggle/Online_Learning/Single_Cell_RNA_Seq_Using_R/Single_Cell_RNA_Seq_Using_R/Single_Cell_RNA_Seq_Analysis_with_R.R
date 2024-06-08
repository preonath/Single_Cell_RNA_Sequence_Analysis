# install.packages("tidyr")
# install.packages("ggplot2")
# install.packages("readr")
#  install.packages("BiocManager")
#BiocManager::install("SingleCellExperiment")

library(tidyr)     # contains tools to tidy data
library(ggplot2) # for plotting
library(readr)     # a package for parsing data




library(SingleCellExperiment)

counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
rownames(counts) <- paste("gene", 1:10, sep = "")
colnames(counts) <- paste("cell", 1:10, sep = "")
sce <- SingleCellExperiment(
  assays = list(counts = counts),
  rowData = data.frame(gene_names = paste("gene_name", 1:10, sep = "")),
  colData = data.frame(cell_names = paste("cell_name", 1:10, sep = ""))
)
sce

normcounts(sce) <- log2(counts(sce) + 1)
sce

dim(normcounts(sce))

head(normcounts(sce))

library(ggplot2)
# library(tidyverse)
library(tidyr)

set.seed(1)
counts <- as.data.frame(matrix(rpois(100, lambda = 10), ncol=10, nrow=10))
Gene_ids <- paste("gene", 1:10, sep = "")
colnames(counts) <- paste("cell", 1:10, sep = "")
counts<-data.frame(Gene_ids, counts)
counts

ggplot(data = counts, mapping = aes(x = cell1, y = cell2))

ggplot(data = counts, mapping = aes(x = cell1, y = cell2)) + geom_point()

counts<-gather(counts, colnames(counts)[2:11], key = 'Cell_ID', value='Counts')
head(counts)

ggplot(counts,aes(x=Cell_ID, y=Counts)) + geom_boxplot()

library(pheatmap)
set.seed(2)
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Cell", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")
pheatmap(test)

library(ggplot2)

if (!require("remotes")) install.packages("remotes")
remotes::install_github('sinhrks/ggfortify')

library(ggfortify)


Principal_Components<-prcomp(test)
autoplot(Principal_Components, label=TRUE)