library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)

# Function to load and preprocess each sample
load_sample <- function(sample_name, project_name, min_features) {
  print(paste("Loading sample:", sample_name))
  
  # Read data from soupX
  sample_data <- Read10X(data.dir = paste0("./matrix_", sample_name, "_soupX"))
  seurat_obj <- CreateSeuratObject(counts = sample_data, project = project_name, min.cells = 3, min.features = min_features)
  
  # Read and apply singlet filtering (doublet removal)
  kept_cells <- read.csv(paste0("./", sample_name, "_kept_cell.csv"), header = TRUE)$x
  seurat_obj_filtered <- subset(seurat_obj, cells = kept_cells)
  
  return(seurat_obj_filtered)
}

# Function to process and save the merged dataset
process_and_save <- function(seurat_obj, output_file) {
  print(paste("Processing and saving:", output_file))
  
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:50)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, min.dist = 0.1, n.neighbors = 30)
  
  saveRDS(seurat_obj, file = output_file)
}

# Define sample groups for merging
sample_groups <- list(
  "proto_bidentis_combined" = list(samples = c("proto_july22_s1", "proto_july22_s2", "proto_apr22_s3"), min_features = 500),
  "nuc_bidentis_combined" = list(samples = c("nuc_fb_s1", "nuc_fb_s2"), min_features = 100),
  "proto_gynandra_combined" = list(samples = c("proto_gy_s1", "proto_gy_s2"), min_features = 200),
  "nuc_gynandra_combined" = list(samples = c("nuc_gy_s1", "nuc_gy_s2"), min_features = 200)
)

# Loop through each group and merge/process samples
for (group_name in names(sample_groups)) {
  samples <- sample_groups[[group_name]]$samples
  min_features <- sample_groups[[group_name]]$min_features
  
  print(paste("Merging samples for:", group_name))
  
  # Load samples dynamically
  sample_list <- lapply(samples, function(sample) load_sample(sample, sample, min_features))
  
  # Merge Seurat objects
  merged_seurat <- Reduce(function(x, y) merge(x, y = y, add.cell.ids = samples), sample_list)
  
  # Apply filtering only for "proto_gynandra_combined"
  if (group_name == "proto_gynandra_combined") {
    print("Applying additional filtering for proto_gynandra_combined")
    merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 500 & nCount_RNA > 1200 & nFeature_RNA < 8000)
  }
  
  # Process and save
  process_and_save(merged_seurat, paste0("./", group_name, ".rds"))
}
