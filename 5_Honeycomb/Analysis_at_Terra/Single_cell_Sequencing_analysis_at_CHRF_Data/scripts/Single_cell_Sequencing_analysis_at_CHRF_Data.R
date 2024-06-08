library(Seurat)
counts <- read.table("/media/asus/275dd380-2319-4638-bcdd-5f65b2b1d4b5/CHRF_Project_Data/Single_Cell/Single_Cell_NP_NextSeq_Batch_24_R/SC_NP_H_7061/4_Quant/SC_NP_H_7061_featurecounts.txt", header = TRUE, row.names = 1)

seuratObject <- CreateSeuratObject(counts = counts)
View(seuratObject)
# Calculate the percent of reads coming from mitochondrial genes
mito.genes <- grep("^MT-", rownames(seuratObject), value = TRUE)
percent.mito <- Matrix::colSums(seuratObject[mito.genes, ]) / Matrix::colSums(seuratObject)
# seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, features = mito.genes)

seuratObject <- subset(seuratObject, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )
seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
# seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
# seuratObject <- ScaleData(seuratObject, vars.to.regress = "percent.mt")
seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
seuratObject <- FindNeighbors(seuratObject, dims = 1:10)
seuratObject <- FindClusters(seuratObject, resolution = 0.5)
seuratObject <- RunUMAP(seuratObject, dims = 1:10)
DimPlot(seuratObject, reduction = "umap")
markers <- FindAllMarkers(seuratObject, only.pos = TRUE, min.pct = 0.25
                          , logfc.threshold = 0.25)

graphql
seuratObject <- RenameIdents(seuratObject, `0` = "CellType1", `1` = "CellType2", ...)
FeaturePlot(seuratObject, features = c("Gene1", "Gene2"))
VlnPlot(seuratObject, features = c("Gene1", "Gene2"))
saveRDS(seuratObject, file = "seuratObject.rds")
write.csv(markers, "markers.csv")