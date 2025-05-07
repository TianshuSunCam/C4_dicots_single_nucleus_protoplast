library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
###bid
bid.combined=readRDS("./bidentis_combined_integrated_05.rds")

###cluster annotation and dimplot
colour_bar=c("#d9d9d9", "#74c476","#238443","#1d91c0",
             "#bf812d","#addd8e", "#f4a582", "#dfc27d",
             "#6baed6", "#9970ab", "#969696","#bdbdbd",
             "#fee391", "#bababa","#d6604d", "#f1b6da",
             "#e0e0e0", "#fee08b",'#fec44f','#238b45','#525252')
DimPlot(bid.combined, reduction = "umap",group.by ='seurat_clusters',label=T,cols=colour_bar)
DefaultAssay(bid.combined)='RNA'
###
#######cluster annotation
PPDK=c('F-bidentis-S0-002372')
PEPC='FBID-00060088' 
NADPME=c('FBID-00003378')
TKL1=c('FBID-00043968')
ACL5=c('FBID-00053265')
VND5=c('F-bidentis-59-000073')
TMO6=c('F-bidentis-4-001300')
OPS=c('F-bidentis-S0-004020')
LACS1=c('FBID-00052318')
KCS10=c('FBID-00008631')
FAMA=c('FBID-00020874')
HTA2=c('F-bidentis-S0-000322')
POK2='FBID-00008020'

query_gene=c(PPDK,PEPC,NADPME,TKL1,
             ACL5,VND5,TMO6,OPS,
             LACS1,KCS10,FAMA,
             HTA2,POK2)

bid.combined$seurat_clusters=factor(bid.combined$seurat_clusters,
                                         levels = c(1,2,5,19,3,8,9,15,4,7,17,12,18,6,14,0,10,11,13,16,20))


DotPlot(object =  bid.combined,
        features = rev(query_gene),scale.by = 'size',group.by = 'seurat_clusters')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14))& scale_colour_gradientn(colours = custom_magma)
#####
