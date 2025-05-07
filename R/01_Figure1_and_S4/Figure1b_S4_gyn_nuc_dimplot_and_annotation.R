library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
nuc_gyn=readRDS(file ='yourpath/nuc_gynandra_combined_ha.rds')
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
Idents(nuc_gyn)='seurat_clusters'
DimPlot(nuc_gyn, reduction = "umap",label=T)
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
LACS1=c('Gyn-10X-HIC-v5-00019253')
PDF2=c('Gyn-10X-HIC-v5-00024581')
FAMA=c('Gyn-10X-HIC-v5-00027865')
HTA2=c('Gyn-10X-HIC-v5-00008303')
POK2='Gyn-10X-HIC-v5-00005859'
query_gene=c(PPDK,PEPC,NAD_ME1,TKL1,
             ACL5,SWEET12,
             PDF2,FAMA,
             HTA2,POK2)
DotPlot(object =  nuc_gyn,
        features = rev(query_gene),scale.by = 'size',group.by = 'seurat_clusters')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14))& scale_colour_gradientn(colours = custom_magma)

###
nuc_gyn=FindSubCluster(nuc_gyn,2 ,subcluster.name = "sub.cluster",resolution = 0.8, graph.name = "RNA_snn")
ACL5=c('Gyn-10X-HIC-v5-00016900')
VND5='Gyn-10X-HIC-v5-00007499'
TMO6='Gyn-10X-HIC-v5-00017445'
OPS='Gyn-10X-HIC-v5-00024240'
SWEET12=c('Gyn-10X-HIC-v5-00016450')
vas_gene=c(ACL5,VND5,TMO6,OPS,SWEET12)
DimPlot(nuc_gyn,group.by = 'sub.cluster',label = T,repel = T)
DotPlot(object =  nuc_gyn,
        features = rev(vas_gene),scale.by = 'size',group.by = 'sub.cluster')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14))& scale_colour_gradientn(colours = custom_magma)
####add celltype
cell_type=nuc_gyn@meta.data$sub.cluster
cell_type=gsub('2_0',"xylem",cell_type)
cell_type=gsub('2_1',"phloem",cell_type)
cell_type=gsub('2_2',"phloem",cell_type)
cell_type=gsub('2_3',"xylem",cell_type)
cell_type=gsub('2_4',"procambium",cell_type)
cell_type=gsub('2_5',"procambium",cell_type)
cell_type=gsub('2_6',"procambium",cell_type)
cell_type=gsub(18,"epidermis",cell_type)
cell_type=gsub(17,"trichome-myb",cell_type)
cell_type=gsub(16,"proliferating",cell_type)
cell_type=gsub(15,"epidermis",cell_type)
cell_type=gsub(14,"xylem",cell_type)
cell_type=gsub(13,"ua",cell_type)
cell_type=gsub(12,"trichome-sty",cell_type)
cell_type=gsub(11,"companion-SE",cell_type)
cell_type=gsub(10,"guard",cell_type)
cell_type=gsub(9,"companion-SE",cell_type)
cell_type=gsub(8,"mesophyll",cell_type)
cell_type=gsub(7,"proliferating",cell_type)
cell_type=gsub(6,"mesophyll",cell_type)
cell_type=gsub(5,"ua",cell_type)
cell_type=gsub(4,"epidermis",cell_type)
cell_type=gsub(3,"epidermis",cell_type)
cell_type=gsub(1,"BS",cell_type)
cell_type=gsub(0,"mesophyll",cell_type)
nuc_gyn <- AddMetaData(object = nuc_gyn, metadata = cell_type, col.name = "cell_type")
nuc_gyn$cell_type=factor(nuc_gyn$cell_type,levels=c('mesophyll','BS',
                                                    'xylem','phloem','procambium','companion-SE',
                                                    'epidermis','guard','trichome-myb','trichome-sty',
                                                    'proliferating','ua'))
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da", "#decbe4", "#de77ae",
             "#dfc27d", "#fec44f","#ffed6f", "#ffed6f",
             "#fbb4ae","#bdbdbd")
DimPlot(nuc_gyn, reduction = "umap",group.by ='cell_type',label=F,cols=colour_bar)

#####
Idents(nuc_gyn)='cell_type'
DotPlot(object =  nuc_gyn,
        features = rev(query_gene),scale.by = 'size',group.by = 'cell_type')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14,angle = 45, vjust = 1, hjust = 1))& scale_colour_gradientn(colours = custom_magma)
#####
DimPlot(nuc_gyn, reduction = "umap",group.by ='orig.ident',label=F)
####
DimPlot(nuc_gyn, reduction = "umap",group.by ='seurat_clusters',label=T)
Idents(nuc_gyn)='cell_type'
nuc_gyn_sub=subset(nuc_gyn, idents=c('mesophyll','BS','xylem','phloem','epidermis','guard','proliferating'))
DimPlot(nuc_gyn_sub, reduction = "umap",group.by ='cell_type',label=T)
query_gene=c(PPDK,PEPC,NAD_ME1,TKL1,
             ACL5,SWEET12,
             PDF2,FAMA,
             HTA2,POK2)
DotPlot(object =  nuc_gyn_sub,
        features = rev(query_gene),scale.by = 'size',group.by = 'cell_type')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14))& scale_colour_gradientn(colours = custom_magma)
