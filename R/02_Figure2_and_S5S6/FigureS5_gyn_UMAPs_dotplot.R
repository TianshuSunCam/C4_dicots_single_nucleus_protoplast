library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
###gyn
gyn.combined=readRDS("./gyn.combined_sub.rds")
DefaultAssay(gyn.combined)='RNA'
###cluster annotation and dimplot
colour_bar=c("#74c476", "#bf812d","#1d91c0","#c994c7",
             "#dfc27d","#238443", "#7fcdbb", "#fe9929",
             "#addd8e", "#8c510a","#f4a582", "#de77ae",
             "#fec44f", "#f768a1","#df65b0", "#d6604d",
             "#f6e8c3", "#fee08b")
DimPlot(gyn.combined, reduction = "umap",label=T,group.by = "seurat_clusters",cols=colour_bar,label.size=6, repel = T)
###
colour_bar=c("#74c476", '#f1b6da',"#bf812d","#fa9fb5",
             "#f4a582", "#de77ae",
             "#fec44f", "#f768a1","#d6604d",
             "#f6e8c3", "#fee08b",
             "#1d91c0","#fcc5c0","#9e9ac8",'#dfc27d','#e7e1ef',
             "#238443","#c994c7", "#7fcdbb", "#8c6bb1",
             "#fe9929",
             "#addd8e", "#8c510a")
DimPlot(gyn.combined, reduction = "umap",label=T,group.by = "new_seurat_clusters",cols=colour_bar,label.size=4, repel = T)
#######cluster annotation
PPDK=c("Gyn-10X-HIC-v5-00028147")
PEPC=c("Gyn-10X-HIC-v5-00018699")
NAD_ME1=c('Gyn-10X-HIC-v5-00027505')
TKL1=c("Gyn-10X-HIC-v5-00019024")
ACL5=c('Gyn-10X-HIC-v5-00016900')
TMO6='Gyn-10X-HIC-v5-00017445'
OPS=c('Gyn-10X-HIC-v5-00024240')
SWEET12=c('Gyn-10X-HIC-v5-00016450')
APL=c('Gyn-10X-HIC-v5-00022559')
LACS1=c('Gyn-10X-HIC-v5-00019253')#Epi
KCS10=c('Gyn-10X-HIC-v5-00028867')
FAMA=c('Gyn-10X-HIC-v5-00027865')
HTA2=c('Gyn-10X-HIC-v5-00008303')
POK2='Gyn-10X-HIC-v5-00005859'
STY1=c('Gyn-10X-HIC-v5-00000133')
Myb23=c('Gyn-10X-HIC-v5-00013351')
query_gene=c(PPDK,PEPC,NAD_ME1,TKL1,
             ACL5,TMO6,OPS,SWEET12,APL,
             LACS1,KCS10,FAMA,
             HTA2,POK2,STY1,Myb23)
gyn.combined$new_seurat_clusters=factor(gyn.combined$new_seurat_clusters,
                                         levels = c("0","5","8","2",
                                                    "5_sub","6_sub","0_sub","2_sub","1_sub","3_sub","4_sub",
                                                    "11","13","1","4","16","7","17","12",
                                                    "10","15","6","9"))

DotPlot(object =  gyn.combined,
        features = rev(query_gene),scale.by = 'size',group.by = 'new_seurat_clusters')+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14))& scale_colour_gradientn(colours = custom_magma)
#####
