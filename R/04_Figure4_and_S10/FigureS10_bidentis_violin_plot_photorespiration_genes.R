library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(scCustomize)
library(ggpubr)
#######
bid.combined.core_nuc=readRDS("./bid.combined.core_nuc.rds")

Idents(bid.combined.core_nuc)="integrated_cell_type"
colour_bar=(c("#238443", "#1d91c0",
              "#9970ab","#f1b6da",
              "#dfc27d", "#fec44f",
              "#fbb4ae"))

#########
DefaultAssay(bid.combined.core_nuc)='RNA'
GLDT='FBID-00047843'
GLDH="F-bidentis-12-000347"
GLDL="FBID-00001892"
SMH1= "F-bidentis-5-000673"
ROG1="F-bidentis-5-001686"
GGAT="F-bidentis-5-001931"
HPR2="F-bidentis-5-000440"
Core_photorespiratory=c(GLDT,GLDH,SMH1,GGAT,HPR2)
Stacked_VlnPlot(seurat_object = bid.combined.core_nuc, features = HPR2,colors_use=colour_bar,x_lab_rotate = TRUE)
DotPlot(  bid.combined.core_nuc, features = rev(Core_photorespiratory),scale.min = 0,
          scale.max = 50)+coord_flip() & scale_colour_gradientn(colours = custom_magma)
