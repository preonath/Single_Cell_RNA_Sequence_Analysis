install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "SummarizedExperiment"),force = TRUE)

ad <- read_h5ad("/home/chrf/Desktop/Data/Single_Cell/Paper_Data/local.h5ad")
