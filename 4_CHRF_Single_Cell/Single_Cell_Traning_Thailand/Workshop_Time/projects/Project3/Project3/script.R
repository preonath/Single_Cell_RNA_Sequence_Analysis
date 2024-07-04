## HCA Africa: Single-Cell RNA-Seq Computational Workshop, March 7, 2023
library(Seurat)
###Project 3 Mouse atlas
setwd("Computational_Projects/Project3")

##Now lets load the mouse data 
mouse.matrix <- read.csv2(file = "counts.csv.gz",sep = ",",row.names = 1)
mouse.matrix[1:5,1:5]
mouse.metadata <- read.csv2(file = "metadata.csv",sep = ",",row.names = 1,header = T)
head(mouse.metadata)
##Lets create the seurat object for Mouse
seuratObject <- CreateSeuratObject(counts = mouse.matrix, meta.data = mouse.metadata, project = "Project3")

##How many cells do we have in each tissues?

## Lets do QC! Use the practical from the previous day and
#check the seurat tutorial! https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

## Once We QC and filter lets do normalization, scaling , PCA, clustering and umap!

## Identify Clusters! Hint: Check for  immunes cells and intestinal cells for example B cells,
#Dendritic Cells, Epithelial cells, Erythrocites, Macrophages,
#Hematopoietic stem progenitor cell, Monocytes , Neutrophils , NK Cells and T Cells
#Use websites like https://singlecell.broadinstitute.org/single_cell and
#https://panglaodb.se/ (if cant connect use an online proxy) as well a protein tissue atlas 

## How many cell do we have in every condition and sample? Is there any batch effect? Do we need to integrate? 
#Tip use the practicals and the seurat tutorial https://satijalab.org/seurat/articles/integration_introduction.html

##How about per cluster/Cell type

##Lets play with Differential exppression

## Lets make pretty plots! 
#Check visualization from Seurat https://satijalab.org/seurat/articles/visualization_vignette.html and GGPlot2
