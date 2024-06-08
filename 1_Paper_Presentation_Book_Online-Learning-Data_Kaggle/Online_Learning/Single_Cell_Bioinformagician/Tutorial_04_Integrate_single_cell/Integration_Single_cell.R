# load libraries
library(Seurat)
library(ggplot2)
#library(tidyverse)
library(gridExtra)

# get data location
dirs <- list.dirs(path = 'data/', recursive = F, full.names = F)

for(x in dirs){
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  cts <- ReadMtx(mtx = paste0('data/',x,'/matrix.mtx.gz'),
                 features = paste0('data/',x,'/features.tsv.gz'),
                 cells = paste0('data/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}
cts

# merge datasets

merged_seurat <- merge(HB17_background, y = c(HB17_PDX, HB17_tumor, HB30_PDX, HB30_tumor, HB53_background,
                                              HB53_tumor),
                       add.cell.ids = ls()[3:9],
                       project = 'HB')


merged_seurat
