library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(scCustomize)
library(ggpubr)

#######
gyn.combined.core_nuc=readRDS("gyn.combined.core_nuc.rds")
Idents(gyn.combined.core_nuc)="new_cell_type"
gyn.combined.core_nuc$new_cell_type=factor(gyn.combined.core_nuc$new_cell_type,
                                       levels=c('mesophyll','BS','xylem','phloem',
                                                'epidermis','guard','proliferating'))
colour_bar=(c("#238443", "#1d91c0",
              "#9970ab","#f1b6da",
              "#dfc27d", "#fec44f",
              "#fbb4ae"))

#########
DefaultAssay(gyn.combined.core_nuc)='RNA'
GLDT= "Gyn-10X-HIC-v5-00007462"
GLDH="Gyn-10X-HIC-v5-00013857"
SHM1= "Gyn-10X-HIC-v5-00020875"
GGAT1="Gyn-10X-HIC-v5-00002976"
SGAT="Gyn-10X-HIC-v5-00003875" 
HPR1="Gyn-10X-HIC-v5-00021074"
Core_photorespiratory=c(GLDT,GLDH,SHM1,GGAT1,HPR1)

custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
DotPlot(  gyn.combined.core_nuc, features = rev(Core_photorespiratory),scale.min = 0,
          scale.max = 50)+coord_flip() & scale_colour_gradientn(colours = custom_magma)
