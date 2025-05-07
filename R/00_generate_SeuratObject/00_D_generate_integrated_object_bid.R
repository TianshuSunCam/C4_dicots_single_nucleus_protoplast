library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)


setwd("youtpath")
nuc_bidentis_combined_ha=readRDS(file = "./nuc_bidentis_combined_ha.rds")
proto_bidentis_combined_ha=readRDS(file = "./proto_bidentis_combined_ha.rds")

bidentis_combined <- merge(x=nuc_bidentis_combined_ha, y = proto_bidentis_combined_ha, 
                           add.cell.ids = c("nuclei", "protoplasts"), project = "bidentis")

rm(nuc_bidentis_combined_ha)
rm(pbmc_proto2_new_filter_single)

bidentis_combined_sample_type=bidentis_combined@meta.data$orig.ident
bidentis_combined_sample_type=gsub('nuc_fb_s1',"nuc",bidentis_combined_sample_type)
bidentis_combined_sample_type=gsub('nuc_fb_s2',"nuc",bidentis_combined_sample_type)
bidentis_combined_sample_type=gsub('proto_july22_s1',"proto",bidentis_combined_sample_type)
bidentis_combined_sample_type=gsub('proto_july22_s2',"proto",bidentis_combined_sample_type)
bidentis_combined_sample_type=gsub('proto_apr22_s3',"proto",bidentis_combined_sample_type)
bidentis_combined= AddMetaData(object = bidentis_combined, metadata = bidentis_combined_sample_type, col.name = "sample_type")

ifnb.list <- SplitObject(bidentis_combined, split.by = "orig.ident")
# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)
# find anchors between samples
bidentis.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# create list of common genes to keep
to_integrate <- Reduce(intersect, lapply(bidentis.anchors@object.list, rownames))
# integrate data and keep full geneset
bidentis.combined <- IntegrateData(anchorset = bidentis.anchors, features.to.integrate = to_integrate)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(bidentis.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
bidentis.combined <- ScaleData(bidentis.combined, verbose = FALSE)
bidentis.combined <- RunPCA(bidentis.combined, npcs = 30, verbose = FALSE)
bidentis.combined <- RunUMAP(bidentis.combined, reduction = "pca", dims = 1:30)
bidentis.combined <- FindNeighbors(bidentis.combined, reduction = "pca", dims = 1:30)
bidentis.combined <- FindClusters(bidentis.combined, resolution = 0.5)

p1=DimPlot(bidentis.combined, reduction = "umap", group.by = "sample_type",cols=c('#fc8d59','#91cf60'))
p2 <- DimPlot(bidentis.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1+p2
DimPlot(bidentis_combined, reduction = "umap",label=T)
saveRDS(bidentis.combined, file = "./bidentis_combined_integrated_05.rds")


