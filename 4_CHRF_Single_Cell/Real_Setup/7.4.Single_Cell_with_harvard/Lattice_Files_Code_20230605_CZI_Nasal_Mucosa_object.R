---
title: "20230605_Create_Object"
output: html_document
date: "2023-06-05"
author: Jaclyn Long
---

Set working directory 
```{r setup, include=FALSE}
knitr::opts_knit$set(root.dir = '~/Pediatric_Swabs')
```

Load necessary libraries 
```{r}
#Essential
library(hdf5r) #Read h5 files 
library(Seurat) #Essential for single cell analysis
 
#Working with dataframes 
#library(readr)
library(dplyr)
library(tibble)

#Plotting 
library(ggplot2)
library(cowplot)
```

Load necessary files 
```{r}
#Define list of pools
pools <- c("20221021_BCH_B003_P002", 
           "20221130_BCH_B005_P002",
           "20221213_BCH_B007_P001",
           "20230222_BCH_B012_P003")


#Load participant metadata table 
metadata <- read.csv(file = "./20230605_lattice_submission_files/20230606_Participant_Metadata_for_scRNAseq_Controls.csv", row.names = "Pool_Fmux_ID")

```

Loop through the list of pools. For each pool, create a seurat object. Then add this to a list of seurat objects. Merge the list of objects into one object. 
Prior to running this chunk, create two folders: 
1. "Alignment_h5_files" and subfolders for each pool containing the filtered_feature_bc_matrix.h5 file.
2. "Freemuxlet", which also has a subfolder for each pool containing the freemuxlet output. The necessary output to generate this object is the file ending in ".clust1.samples.gz". 
```{r}
#create empty list
obj.list <- list()

for(pool in pools){
  print(pool)
  #read h5 file 
  h5.path <- paste("./Alignment_h5_files/", pool, "/filtered_feature_bc_matrix.h5", sep = "")
  obj.pool <- Read10X_h5(h5.path)
  
  #create seurat object 
  obj.pool <- CreateSeuratObject(obj.pool)
  
  obj.pool$Pool <- pool
  
  #print seurat object info 
  print(obj.pool)
  
  #load freemuxlet results 
  fmux.path <- paste("./Freemuxlet/", pool, "/", pool, ".clust1.samples.gz", sep = "")
  fmux <- as.data.frame(read.table(gzfile(fmux.path), header = TRUE))
  
  #format table, subset to only cells in the seurat object and sort by cell barcode 
  rownames(fmux) <- fmux$BARCODE
  fmux <- subset(fmux, subset = fmux$BARCODE %in% colnames(obj.pool))
  fmux <- fmux[order(fmux$BARCODE),]
  
  #add freemuxlet results to seurat object as metadata 
  obj.pool <- subset(obj.pool, cells = rownames(fmux))
  obj.pool <- AddMetaData(obj.pool, metadata = fmux$DROPLET.TYPE, col.name = "Fmux_Droplet_Type")
  obj.pool <- AddMetaData(obj.pool, metadata = fmux$SNG.BEST.GUESS, col.name = "Fmux_Cluster_ID")
  
  #add metadata column with site, batch, pool, and freemuxlet ID 
  obj.pool <- AddMetaData(obj.pool, paste(obj.pool$Pool, obj.pool$Fmux_Cluster_ID, sep = "_"), "Pool_FmuxID") 
  
  #Filter to droplets defined as singlets by Freemuxlet 
  obj.pool <- subset(obj.pool, Fmux_Droplet_Type == "SNG")
  
  #split object by fmux ID
  obj.pool.list <- SplitObject(obj.pool, split.by = "Pool_FmuxID")

  #Add metatadata from participant table
  for (obj.pt in obj.pool.list){
    pt <- unique(obj.pt$Pool_FmuxID)
    for (var in colnames(metadata)){
      obj.pt <- AddMetaData(obj.pt, metadata[pt, var], col.name = var)
    }

    obj.pool.list[[pt]] <- obj.pt
  }
  
  #merge split objects back into one object for the pool 
  obj.pool <- merge(x = obj.pool.list[[1]], y = obj.pool.list[2:length(obj.pool.list)], add.cell.ids = names(obj.pool.list))
  
  #print number of cells from each participant 
  print(table(obj.pool$Participant_ID))
  
  #remove unnecessary metadata
  obj.pool$Pool <- NULL
  obj.pool$Pool_FmuxID <- NULL 
  obj.pool$Pool_Participant_ID <- NULL 
  
  #add seurat object to list 
  obj.list[[pool]] <- obj.pool
  
}

#Merge list into one seurat object
obj <- merge(x = obj.list[[1]], y = obj.list[2:length(obj.list)], add.cell.ids = pools)
obj

rm(obj.list, obj.pool, obj.pool.list, obj.pt)
```

QC plots 
```{r}
#add percent mitochondrial
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

#plot by pool
ggplot(as.data.frame(table(obj$Pool_ID)) %>% rename(Pool_ID = Var1, Number_Cells = Freq), 
       aes(x = Pool_ID, y = Number_Cells, fill = Pool_ID)) +
  geom_bar(stat = "identity", width = .6, color = "black") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,20000)) + 
  ylab("# of Cells") +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 90, vjust = .9), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14), 
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  geom_text(aes(label = Number_Cells), vjust = -0.2, size = 3)

VlnPlot(obj, features = "nFeature_RNA", log = TRUE, group.by = "Pool_ID") * 
  theme(axis.title.x = element_blank(), legend.position = "none") + 
  geom_hline(yintercept = 300, color = "red") 

VlnPlot(obj, features = "nCount_RNA", log = TRUE, group.by = "Pool_ID") * 
  theme(axis.title.x = element_blank(), legend.position = "none") + 
  geom_hline(yintercept = 500, color = "red")

VlnPlot(obj, features = "percent.mt",  group.by = "Pool_ID") * 
  theme(axis.title.x = element_blank(), legend.position = "none") + 
  geom_hline(yintercept = 35, color = "red")



#plot by participant

ggplot(as.data.frame(table(obj$Sample_ID)) %>% rename(Sample_ID = Var1, Number_Cells = Freq), 
       aes(x = Sample_ID, y = Number_Cells, fill = Sample_ID)) +
  geom_bar(stat = "identity", width = .6, color = "black") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,8500)) + 
  ylab("# of Cells") +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 90, vjust = .9), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14), 
        plot.margin = unit(c(.5,.5,.5,.5), "cm")) +
  geom_text(aes(label = Number_Cells), vjust = -0.2, size = 3)


VlnPlot(obj, features = "nFeature_RNA", log = TRUE, group.by = "Sample_ID") + 
  theme(axis.title.x = element_blank(), legend.position = "none") + 
  geom_hline(yintercept = 300, color = "red") 

VlnPlot(obj, features = "nCount_RNA", log = TRUE, group.by = "Sample_ID") +
  theme(axis.title.x = element_blank(), legend.position = "none") + 
  geom_hline(yintercept = 500, color = "red")

VlnPlot(obj, features = "percent.mt",  group.by = "Sample_ID") + 
  theme(axis.title.x = element_blank(), legend.position = "none") + 
  geom_hline(yintercept = 35, color = "red")

```

Preprocessing(filtering, normalization, PCA)
```{r}
#Filter based on QC metrics 
obj <- subset(obj, subset = nFeature_RNA > 300 & nCount_RNA > 500 & percent.mt < 35)
table(obj$Sample_ID)

#Run SCTransform normalization
obj <- SCTransform(obj, vars.to.regress = "percent.mt")

#PCA 
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
ElbowPlot(obj, ndims = 40)
```

Run clustering and generate UMAP 
```{r}

#clustering 
obj <- FindNeighbors(obj, dims = 1:26) 
obj <- FindClusters(obj, resolution = 0.4) 
obj <- RunUMAP(obj, dims = 1:26)
DimPlot(obj, label = TRUE)
DimPlot(obj, group.by = "Pool_ID", shuffle = TRUE)
DimPlot(obj, group.by = "Sample_ID", shuffle = TRUE)

```

Add coarse cell annotations 
```{r}
obj <- AddMetaData(obj, NA, col.name = "Coarse_Cell_Type")

obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 0 ] <- "Secretory/Goblet" 
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 1 ] <- "Ciliated" 
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 2 ] <- "Ciliated"  
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 3 ] <- "Basal/Club" 
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 4 ] <- "T/NK" 
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 5 ] <- "Secretory/Goblet"
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 6 ] <- "Ciliated" 
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 7 ] <- "Ciliated"  
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 8 ] <- "Myeloid/APC"  
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 9 ] <- "Myeloid/APC"  
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 10 ] <- "Secretory/Goblet"
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 11 ] <- "Ionocytes" 
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 12 ] <- "T/NK"
obj@meta.data$Coarse_Cell_Type[obj@meta.data$seurat_clusters == 13 ] <- "Squamous" 
 

cluster.order <- c(
  "Basal/Club",
  "Secretory/Goblet",
  "Squamous",
  "Ionocytes",
  "Ciliated",
  "T/NK",
  "Myeloid/APC" 
)

obj@meta.data$Coarse_Cell_Type<- factor(obj@meta.data$Coarse_Cell_Type, levels = cluster.order)



DimPlot(obj, group.by = "Coarse_Cell_Type", shuffle = TRUE) + 
  labs(title = "Cell Types") + 
  theme(axis.text = element_blank(), axis.ticks= element_blank())

```

Stacked violin Plot of cell type marker genes 
```{r}
marker.genes <- c(
 #basal
"KRT5",
"KRT15", 

#club 
"SCGB1A1",
"BPIFA1",

#goblet 
"MUC5AC",
"MUC2", 
"AQP5",

#squamous
"SPRR1B",
"SPRR2A",
"KRT17", 
"KRT23", 

#ionocytes
"RARRES2",
"SCNN1B",
"CFTR",

#ciliated 
"CAPS",
"DNAAF1",
"DNAH12", 
"TPPP3",

#T
"CD3D",
"TRAC", 
"IFNG",
"NKG7", 

#myeloid 
"TYROBP",
"LYZ",

#HLA
"HLA-DPA1",
"HLA-DPB1",
"HLA-DQA1"

)


VlnPlot(obj, features = marker.genes, group.by = "Coarse_Cell_Type", fill.by = "ident",  stack = TRUE, pt.size = 0) + 
  theme(legend.position = "none", axis.title.y = element_blank()) + 
  labs(title = "Cell Type Marker Genes")
```

Save object 
```{r}
saveRDS(obj, file = "./20230605_lattice_submission_files/20230606_CZI_Nasal_Mucosa.rds")

obj <- readRDS("./20230605_lattice_submission_files/20230606_CZI_Nasal_Mucosa.rds")

```

