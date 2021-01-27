library(Seurat)
library(tidyverse)
library(SCINA)

Donor1 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138939')
Donor2 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138940')
Donor3 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138941')
Donor4 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138942')
Donor5 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138943')
Donor6 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138944')
Donor7 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138945')
Donor8 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138946')
Donor9 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138947')
Donor10 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138948')
Donor11 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138949')
Donor12 <- Read10X('D:/Endothelial scRNA-Seq project/PMID 30380031 pancreatic islets/GSM3138950')

#Create Seurat objects
Donor1 <- CreateSeuratObject(counts = Donor1, project = "Donor1", min.cells = 3, min.features = 200)
Donor2 <- CreateSeuratObject(counts = Donor2, project = "Donor2", min.cells = 3, min.features = 200)
Donor3 <- CreateSeuratObject(counts = Donor3, project = "Donor3", min.cells = 3, min.features = 200)
Donor4 <- CreateSeuratObject(counts = Donor4, project = "Donor4", min.cells = 3, min.features = 200)
Donor5 <- CreateSeuratObject(counts = Donor5, project = "Donor5", min.cells = 3, min.features = 200)
Donor6 <- CreateSeuratObject(counts = Donor6, project = "Donor6", min.cells = 3, min.features = 200)
Donor7 <- CreateSeuratObject(counts = Donor7, project = "Donor7", min.cells = 3, min.features = 200)
Donor8 <- CreateSeuratObject(counts = Donor8, project = "Donor8", min.cells = 3, min.features = 200)
Donor9 <- CreateSeuratObject(counts = Donor9, project = "Donor9", min.cells = 3, min.features = 200)
Donor10 <- CreateSeuratObject(counts = Donor10, project = "Donor10", min.cells = 3, min.features = 200)
Donor11 <- CreateSeuratObject(counts = Donor11, project = "Donor11", min.cells = 3, min.features = 200)
Donor12 <- CreateSeuratObject(counts = Donor12, project = "Donor12", min.cells = 3, min.features = 200)

#Create list of samples
Pancreas.list <- list(Donor1, Donor2, Donor3, Donor4, Donor5, Donor6, Donor7, Donor8, Donor9, Donor10, Donor11, Donor12)

#Compute mito %
for(i in 1:length(Pancreas.list)) {
  Pancreas.list[[i]][["percent.mito"]] <- PercentageFeatureSet(Pancreas.list[[i]], pattern = "^MT-")
}


#Subset
for(i in 1:length(Pancreas.list)) {
  Pancreas.list[[i]] <- subset(Pancreas.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mito < 30)
}

#Merge datasets
Pancreas <- merge(Pancreas.list[[1]], y = c(Pancreas.list[[2]],Pancreas.list[[3]],Pancreas.list[[4]],Pancreas.list[[5]],Pancreas.list[[6]],Pancreas.list[[7]],Pancreas.list[[8]],Pancreas.list[[9]],Pancreas.list[[10]],Pancreas.list[[11]],Pancreas.list[[12]]), add.cell.ids = c("Donor1", "Donor2", "Donor3", "Donor4", "Donor5", "Donor6", "Donor7", "Donor8", "Donor9", "Donor10", "Donor11", "Donor12"))

#Standard Seurat workflow
Pancreas <- NormalizeData(Pancreas, normalization.method = "LogNormalize", scale.factor = 10000)
Pancreas <- FindVariableFeatures(Pancreas, selection.method = "vst", nfeatures = 2000)
Pancreas <- ScaleData(Pancreas)
Pancreas <- RunPCA(Pancreas, features = VariableFeatures(object = Pancreas))
Pancreas <- FindNeighbors(Pancreas, dims = 1:10)
Pancreas <- FindClusters(Pancreas, resolution = 0.5)
Pancreas <- RunUMAP(Pancreas, dims = 1:10)

DimPlot(Pancreas, reduction = "umap", label = TRUE)
VlnPlot(Pancreas, features = c('CLDN5', 'VWF', 'FLT1', 'PECAM1'), ncol = 2)
Pancreas <- subset(Pancreas, idents = c('10'))

#Set list of endothelial markers
Endothelial_Markers <- list(c('CLDN5', 'VWF', 'FLT1', 'PECAM1'))
names(Endothelial_Markers) <- 'Endothelial cells'

#Run SCINA for cell type prediction
SCINA_results <- SCINA(Pancreas@assays$RNA@data,
                       Endothelial_Markers,
                       max_iter = 2000, 
                       convergence_n = 100, 
                       convergence_rate = 0.999, 
                       sensitivity_cutoff = 0.9, 
                       rm_overlap=FALSE, 
                       allow_unknown=TRUE)
Pancreas$cell_labels <- SCINA_results$cell_labels
DimPlot(Pancreas,reduction = "umap", pt.size = 1, label = TRUE, group.by = 'cell_labels')

#Subset for endothelial and write object
Pancreas <- subset(Pancreas, cell_labels == 'Endothelial cells')
Pancreas$organ <- 'Pancreas'
write_rds(Pancreas, 'Pancreas.rds')
