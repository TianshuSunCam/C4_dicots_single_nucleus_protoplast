library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)

proto_gyn=readRDS(file ='yourpath/proto_gynandra_combined_ha.rds')
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
Idents(proto_gyn)='seurat_clusters'
DimPlot(proto_gyn, reduction = "umap",label=T)

#######cluster annotation
PPDK=c("Gyn-10X-HIC-v5-00028147")
PEPC=c("Gyn-10X-HIC-v5-00018699")
NAD_ME1=c('Gyn-10X-HIC-v5-00027505')
TKL1=c("Gyn-10X-HIC-v5-00019024")
ACL5=c('Gyn-10X-HIC-v5-00016900')
TMO6='Gyn-10X-HIC-v5-00017445'
SWEET12=c('Gyn-10X-HIC-v5-00016450')
APL=c('Gyn-10X-HIC-v5-00022559')
KCS10=c('Gyn-10X-HIC-v5-00028867')
FAMA=c('Gyn-10X-HIC-v5-00027865')
HTA2=c('Gyn-10X-HIC-v5-00008303')
POK2='Gyn-10X-HIC-v5-00005859'
query_gene=c(PPDK,PEPC,NAD_ME1,TKL1,
             ACL5,TMO6,
             KCS10,FAMA,
             HTA2,POK2)
DotPlot(object =  proto_gyn,
        features = rev(query_gene),scale.by = 'size',group.by = 'seurat_clusters')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14))& scale_colour_gradientn(colours = custom_magma)
###
proto_gyn=FindSubCluster(proto_gyn,10 ,subcluster.name = "sub.cluster",resolution = 0.4, graph.name = "RNA_snn")
ACL5=c('Gyn-10X-HIC-v5-00016900')
VND5='Gyn-10X-HIC-v5-00007499'
TMO6='Gyn-10X-HIC-v5-00017445'
OPS='Gyn-10X-HIC-v5-00024240'
SWEET12=c('Gyn-10X-HIC-v5-00016450')
vas_gene=c(ACL5,VND5,TMO6,OPS,SWEET12)
DimPlot(proto_gyn,group.by = 'sub.cluster',label = T,repel = T)
DotPlot(object =  proto_gyn,
        features = rev(vas_gene),scale.by = 'size',group.by = 'sub.cluster')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14))& scale_colour_gradientn(colours = custom_magma)
####add celltype
cell_type=proto_gyn@meta.data$sub.cluster
cell_type=gsub('10_0',"xylem",cell_type)
cell_type=gsub('10_1',"phloem",cell_type)
cell_type=gsub('10_2',"xylem",cell_type)
cell_type=gsub(14,"ua",cell_type)
cell_type=gsub(13,"companion-SE",cell_type)
cell_type=gsub(12,"epidermis",cell_type)
cell_type=gsub(11,"proliferating",cell_type)
cell_type=gsub(9,"epidermis",cell_type)
cell_type=gsub(8,"BS",cell_type)
cell_type=gsub(7,"ua",cell_type)
cell_type=gsub(6,"epidermis",cell_type)
cell_type=gsub(5,"mesophyll",cell_type)
cell_type=gsub(4,"mesophyll",cell_type)
cell_type=gsub(3,"guard",cell_type)
cell_type=gsub(2,"epidermis",cell_type)
cell_type=gsub(1,"epidermis",cell_type)
cell_type=gsub(0,"mesophyll",cell_type)
proto_gyn <- AddMetaData(object = proto_gyn, metadata = cell_type, col.name = "cell_type")
proto_gyn$cell_type=factor(proto_gyn$cell_type,levels=c('mesophyll','BS',
                                                    'xylem','phloem','companion-SE',
                                                    'epidermis','guard',
                                                    'proliferating','ua'))
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da", "#de77ae",
             "#dfc27d", "#fec44f",
             "#fbb4ae","#bdbdbd")

DimPlot(proto_gyn, reduction = "umap",group.by ='cell_type',label=F,cols=colour_bar)
#####
Idents(proto_gyn)='cell_type'
DotPlot(object =  proto_gyn,
        features = rev(query_gene),scale.by = 'size',group.by = 'cell_type')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1))& scale_colour_gradientn(colours = custom_magma)
####
DimPlot(proto_gyn, reduction = "umap",group.by ='orig.ident',label=F)
####
DimPlot(proto_gyn, reduction = "umap",group.by ='seurat_clusters',label=F)
Idents(proto_gyn)='cell_type'
proto_gyn_sub=subset(proto_gyn, idents=c('mesophyll','BS','xylem','phloem','epidermis','guard','proliferating'))
DimPlot(proto_gyn_sub, reduction = "umap",group.by ='cell_type',label=T)
query_gene=c(PPDK,PEPC,NAD_ME1,TKL1,
             ACL5,SWEET12,
             PDF2,FAMA,
             HTA2,POK2)
DotPlot(object =  proto_gyn_sub,
        features = rev(query_gene),scale.by = 'size',group.by = 'cell_type')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14))& scale_colour_gradientn(colours = custom_magma)
