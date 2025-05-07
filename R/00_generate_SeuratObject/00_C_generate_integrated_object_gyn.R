library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)

setwd("yourpath")
nuc_gynandra_combined_ha=readRDS("./nuc_gynandra_combined_ha.rds")
proto_gynandra_combined_ha=readRDS("./proto_gynandra_combined_ha.rds")
gyn_combined <- merge(nuc_gynandra_combined_ha, y = proto_gynandra_combined_ha, 
                      add.cell.ids = c("nuclei", "protoplasts"), project = "Gynandra")
gyn_combined_sample_type=gyn_combined@meta.data$orig.ident
gyn_combined_sample_type=gsub('nuc_gy_s1',"nuc",gyn_combined_sample_type)
gyn_combined_sample_type=gsub('nuc_gy_s2',"nuc",gyn_combined_sample_type)
gyn_combined_sample_type=gsub('proto_gy_s1',"proto",gyn_combined_sample_type)
gyn_combined_sample_type=gsub('proto_gy_s2',"proto",gyn_combined_sample_type)
gyn_combined= AddMetaData(object = gyn_combined, metadata = gyn_combined_sample_type, col.name = "sample_type")

ifnb.list <- SplitObject(gyn_combined, split.by = "sample_type")
# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)
gyn.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
gyn.combined <- IntegrateData(anchorset = gyn.anchors)
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(gyn.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
gyn.combined <- ScaleData(gyn.combined, verbose = FALSE)
gyn.combined <- RunPCA(gyn.combined, npcs = 30, verbose = FALSE)
gyn.combined <- RunUMAP(gyn.combined, reduction = "pca", dims = 1:30)
gyn.combined <- FindNeighbors(gyn.combined, reduction = "pca", dims = 1:30)
gyn.combined <- FindClusters(gyn.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(gyn.combined, reduction = "umap", group.by = "sample_type")
p2 <- DimPlot(gyn.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
colour_bar=c("#ccebc5", "#b3cde3", "#decbe4",
             "#fddaec", "#e5d8bd", "#fed9a6",
             "#ffed6f", "#ffed6f","#fbb4ae", "#f2f2f2")
DimPlot(gyn.combined, reduction = "umap",split.by ='sample_type',group.by ='cell_type_id',label=T,cols=colour_bar)
DimPlot(gyn_combined, reduction = "umap",label=T)

integrated_cell_type=gyn.combined@meta.data$seurat_clusters
integrated_cell_type=gsub(17,"trichome-myb",integrated_cell_type)
integrated_cell_type=gsub(16,"epidermis",integrated_cell_type)
integrated_cell_type=gsub(15,"proliferating",integrated_cell_type)
integrated_cell_type=gsub(14,"vascular",integrated_cell_type)
integrated_cell_type=gsub(13,"companion-SE",integrated_cell_type)
integrated_cell_type=gsub(12,"trichome-sty",integrated_cell_type)
integrated_cell_type=gsub(11,"companion-SE",integrated_cell_type)
integrated_cell_type=gsub(10,"proliferating",integrated_cell_type)

###cluster >=10 change first
integrated_cell_type=gsub(9,"ua",integrated_cell_type)
integrated_cell_type=gsub(8,"mesophyll",integrated_cell_type)
integrated_cell_type=gsub(7,"guard",integrated_cell_type)
integrated_cell_type=gsub(6,"ua",integrated_cell_type)
integrated_cell_type=gsub(5,"mesophyll",integrated_cell_type)
integrated_cell_type=gsub(4,"epidermis",integrated_cell_type)
integrated_cell_type=gsub(3,"vascular",integrated_cell_type)
integrated_cell_type=gsub(2,"BS",integrated_cell_type)
#####cluster 0 and cluster 1 should change last
integrated_cell_type=gsub(1,"epidermis",integrated_cell_type)
integrated_cell_type=gsub(0,"mesophyll",integrated_cell_type)

gyn.combined <- AddMetaData(object = gyn.combined, 
                            metadata = integrated_cell_type, col.name = "integrated_cell_type")
gyn.combined$integrated_cell_type <- factor(gyn.combined$integrated_cell_type,
                                            levels=c('mesophyll','BS','vascular',
                                                     'companion-SE','epidermis','guard',
                                                     'trichome-sty','trichome-myb',
                                                     'proliferating','ua'))

gyn.combined$integrated_CT_marker<- factor(gyn.combined$integrated_CT_marker,levels=c('mesophyll','BS','vascular',
                                                                                      'companion','SE','epidermis','guard',
                                                                                      'trichome-sty','trichome-myb',
                                                                                      'proliferating','ua-6','ua-9'))
DimPlot(gyn.combined, reduction = "umap",group.by ='integrated_CT_marker',label=T,cols=colour_bar)

Idents(gyn.combined)='integrated_cell_type'
vasculature =subset(x = gyn.combined, idents = c("vascular"))
head(vasculature)
DimPlot(vasculature, reduction = "umap",group.by ='seurat_clusters',label=T)
vasculature_reclustering=vasculature
vasculature_reclustering <- FindNeighbors(object = vasculature_reclustering, dims = 1:30)
vasculature_reclustering <- FindClusters(object = vasculature_reclustering, resolution = 0.5)
vasculature_reclustering <- RunUMAP(vasculature_reclustering, dims = 1:10,min.dist = 0.1,n.neighbors=30)
DimPlot(vasculature_reclustering, reduction = "umap",label=T,group.by = 'seurat_clusters',label.size = 6)
######
integrated_subtype=vasculature_reclustering@meta.data$seurat_clusters
integrated_subtype=gsub(6,"xylem",integrated_subtype)
integrated_subtype=gsub(5,"xylem",integrated_subtype)
integrated_subtype=gsub(4,"procambium",integrated_subtype)
integrated_subtype=gsub(3,"procambium",integrated_subtype)
integrated_subtype=gsub(2,"phloem",integrated_subtype)
integrated_subtype=gsub(1,"procambium",integrated_subtype)
integrated_subtype=gsub(0,"phloem",integrated_subtype)
vasculature_reclustering <- AddMetaData(object = vasculature_reclustering, 
                                        metadata = integrated_subtype, col.name = "sub_cell_type")
DimPlot(vasculature_reclustering, reduction = "umap",label=T,group.by = 'sub_cell_type',label.size = 6)

#####
integrated_subcluster=vasculature_reclustering@meta.data$seurat_clusters
integrated_subcluster=gsub(6,"6_sub",integrated_subcluster)
integrated_subcluster=gsub(5,"5_sub",integrated_subcluster)
integrated_subcluster=gsub(4,"4_sub",integrated_subcluster)
integrated_subcluster=gsub(3,"3_sub",integrated_subcluster)
integrated_subcluster=gsub(2,"2_sub",integrated_subcluster)
integrated_subcluster=gsub(1,"1_sub",integrated_subcluster)
integrated_subcluster=gsub(0,"0_sub",integrated_subcluster)
vasculature_reclustering <- AddMetaData(object = vasculature_reclustering, 
                                        metadata = integrated_subcluster, col.name = "sub_seurat_cluster")
DimPlot(vasculature_reclustering, reduction = "umap",label=T,group.by = 'sub_seurat_cluster',label.size = 6)
gyn.combined_sub=gyn.combined
Cells <- WhichCells(gyn.combined_sub)
head(vasculature_reclustering)
new_cell_type=character()
for(i in 1:length(Cells)){
  if (Cells[i] %in% rownames(vasculature_reclustering@meta.data) ){
    new_cell_type[i]=as.character(vasculature_reclustering@meta.data[which(rownames(vasculature_reclustering@meta.data)==Cells[i]),]$sub_cell_type)
  } else {
    new_cell_type[i]=as.character(gyn.combined_sub@meta.data[which(rownames(gyn.combined_sub@meta.data)==Cells[i]),]$integrated_cell_type)
  } 
}

gyn.combined_sub <- AddMetaData(object = gyn.combined_sub, metadata = new_cell_type, col.name = "new_cell_type")
colorBar=c('#F8766D','#B79F00','#00BA38','#00BFC4',
           '#619CFF','#F564E3','#cccccc')
DimPlot(gyn.combined_sub, reduction = "umap",label=T,group.by = "new_cell_type")

#######
Cells <- WhichCells(gyn.combined_sub)
head(vasculature_reclustering)
new_cell_cluster=character()
for(i in 1:length(Cells)){
  if (Cells[i] %in% rownames(vasculature_reclustering@meta.data) ){
    new_cell_cluster[i]=as.character(vasculature_reclustering@meta.data[which(rownames(vasculature_reclustering@meta.data)==Cells[i]),]$sub_seurat_cluster)
  } else {
    new_cell_cluster[i]=as.character(gyn.combined_sub@meta.data[which(rownames(gyn.combined_sub@meta.data)==Cells[i]),]$seurat_clusters)
  } 
}

gyn.combined_sub <- AddMetaData(object = gyn.combined_sub, metadata = new_cell_cluster, col.name = "new_seurat_clusters")

DimPlot(gyn.combined_sub, reduction = "umap",label=T,group.by = "new_seurat_clusters")

gyn.combined_sub$new_cell_type<- factor(gyn.combined_sub$new_cell_type,
                                        levels=c('mesophyll','BS',
                                                 'xylem','phloem','procambium','companion-SE',
                                                 'epidermis','guard','trichome-myb','trichome-sty',
                                                 'proliferating','ua'))
saveRDS(gyn.combined_sub, file = "./gyn.combined_sub.rds")