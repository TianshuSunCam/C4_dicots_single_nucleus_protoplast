library(dplyr)
library(SoupX)
library(Seurat)
library(DropletUtils)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(viridis)

# Load sample table
sample_table <- read.csv("doublet_table.csv", sep = "\t", header = TRUE)

for (i in 1:nrow(sample_table)) {
  project_name <- sample_table$project_name[i]
  DoubletRate <- sample_table$DoubletRate[i]
  
  print(paste("Processing:", project_name))
  
  
  # Load data
  sc <- load10X(paste0("path/", project_name, "/outs"))
  sc <- autoEstCont(sc)
  out <- adjustCounts(sc)
  output_dir <- paste0("./matrix_", project_name, "_soupX")
  DropletUtils:::write10xCounts(output_dir, out, version = "3")
  
  # Setup Seurat object
  sample_data <- Read10X(data.dir = output_dir)
  souped_sample <- CreateSeuratObject(counts = sample_data, project = project_name, min.cells = 3, min.features = 200)
  
  # Normalize data
  souped_sample <- NormalizeData(souped_sample, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identify highly variable features
  souped_sample <- FindVariableFeatures(souped_sample, selection.method = "vst", nfeatures = 2000)
  
  # Scale data
  all.genes <- rownames(souped_sample)
  souped_sample <- ScaleData(souped_sample, features = all.genes)
  
  # PCA
  souped_sample <- RunPCA(souped_sample, features = VariableFeatures(object = souped_sample))
  
  # Clustering
  souped_sample <- FindNeighbors(souped_sample, dims = 1:30)
  souped_sample <- FindClusters(souped_sample, resolution = 0.5)
  
  # UMAP
  souped_sample <- RunUMAP(souped_sample, dims = 1:10, min.dist = 0.1, n.neighbors = 30)
  
  # Save full dataset
  saveRDS(souped_sample, file = paste0("./", project_name, "_souped.rds"))
  
  # Doublet Finder
  pc.num <- 1:10
  sweep.res.list <- paramSweep_v3(souped_sample, PCs = pc.num, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  
  homotypic.prop <- modelHomotypic(souped_sample$seurat_clusters)
  nExp_poi <- round(DoubletRate * ncol(souped_sample))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Finding doublets
  souped_sample <- doubletFinder_v3(souped_sample, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                                    nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
  
  DimPlot(souped_sample, group.by = colnames(souped_sample@meta.data)[7], cols = c("black", "gold", "red"), reduction = "umap")
  
  # Extract singlets
  Idents(souped_sample) <- colnames(souped_sample@meta.data)[7]
  souped_sample_single <- subset(souped_sample, 
                                 cells = rownames(souped_sample@meta.data)[which(souped_sample@meta.data[,7] == 'Singlet')])
  
  # Save singlets
  saveRDS(souped_sample_single, file = paste0("./", project_name, "_singlet.rds"))
  
  # Save singlet cell names
  sample_kept_cell <- rownames(souped_sample@meta.data)[which(souped_sample@meta.data[,7] == 'Singlet')]
  write.csv(x = sample_kept_cell, file = paste0("./", project_name, "_kept_cell.csv"), row.names = FALSE, quote = FALSE)
}

print("Processing complete for all samples!")