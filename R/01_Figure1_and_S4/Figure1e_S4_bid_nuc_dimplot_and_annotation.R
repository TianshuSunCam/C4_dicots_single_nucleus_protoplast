library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)

nuc_bid=readRDS(file ='yourpath/nuc_bidentis_combined_ha.rds')
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
Idents(nuc_bid)='seurat_clusters'
DimPlot(nuc_bid, reduction = "umap",label=T)

#######cluster annotation
PPDK=c('F-bidentis-S0-002372')
PEPC='FBID-00060088' #'Fbid_CUFF.4218'
NADPME=c('FBID-00003378')
TKL1=c('FBID-00043968')
ACL5=c('FBID-00053265')
TMO6=c('F-bidentis-4-001300')
OPS=c('F-bidentis-S0-004020','F-bidentis-S0-001349')# Query: AT3G09070 'Gyn-10X-HIC-v5-00011689',
KCS10=c('FBID-00008631')
FAMA=c('FBID-00020874')
HTA2=c('F-bidentis-S0-000322')
POK2='FBID-00008020'#AT3G19050.1
query_gene=c(PPDK,PEPC,NADPME,TKL1,
             ACL5,TMO6,
             KCS10,FAMA,
             HTA2,POK2)
DotPlot(object =  nuc_bid,
        features = rev(query_gene),scale.by = 'size',group.by = 'seurat_clusters')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14))& scale_colour_gradientn(colours = custom_magma)
####add celltype
cell_type=nuc_bid@meta.data$seurat_clusters
cell_type=gsub(18,"ua",cell_type)
cell_type=gsub(17,"guard",cell_type)
cell_type=gsub(16,"proliferating",cell_type)
cell_type=gsub(15,"ua",cell_type)
cell_type=gsub(14,"ua",cell_type)
cell_type=gsub(13,"mesophyll",cell_type)
cell_type=gsub(12,"epidermis",cell_type)
cell_type=gsub(11,"proliferating",cell_type)
cell_type=gsub(10,"proliferating",cell_type)
cell_type=gsub(9,"ua",cell_type)
cell_type=gsub(8,"ua",cell_type)
cell_type=gsub(7,"epidermis",cell_type)
cell_type=gsub(6,"phloem",cell_type)
cell_type=gsub(5,"xylem",cell_type)
cell_type=gsub(4,"proliferating",cell_type)
cell_type=gsub(3,"BS",cell_type)
cell_type=gsub(2,"epidermis",cell_type)
cell_type=gsub(1,"mesophyll",cell_type)
cell_type=gsub(0,"ua",cell_type)
nuc_bid <- AddMetaData(object = nuc_bid, metadata = cell_type, col.name = "cell_type")
nuc_bid$cell_type=factor(nuc_bid$cell_type,levels=c(c('mesophyll','BS','xylem',"phloem",
                                                      'epidermis','guard',
                                                      'proliferating','ua')))
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da",
             "#dfc27d", "#fec44f",
             "#fbb4ae","#bdbdbd")
DimPlot(nuc_bid, reduction = "umap",group.by ='cell_type',label=F,cols=colour_bar)
#####
Idents(nuc_bid)='cell_type'
DotPlot(object =  nuc_bid,
        features = rev(query_gene),scale.by = 'size',group.by = 'cell_type')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14))& scale_colour_gradientn(colours = custom_magma)
#####
DimPlot(nuc_bid, reduction = "umap",group.by ='orig.ident',label=F)
DimPlot(nuc_bid, reduction = "umap",group.by ='seurat_clusters',label=T)
Idents(nuc_bid)='cell_type'
nuc_bid_sub=subset(nuc_bid, idents=c('mesophyll','BS','xylem','phloem','epidermis','guard','proliferating'))
DimPlot(nuc_bid_sub, reduction = "umap",group.by ='cell_type',label=T)
query_gene=c(PPDK,PEPC,NADPME,TKL1,
             ACL5,TMO6,
             KCS10,FAMA,
             HTA2,POK2)
DotPlot(object =  nuc_bid_sub,
        features = rev(query_gene),scale.by = 'size',group.by = 'cell_type')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14))& scale_colour_gradientn(colours = custom_magma)
