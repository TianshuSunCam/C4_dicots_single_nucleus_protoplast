library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(scCustomize)
library(ggpubr)
#######
gyn.combined_sub=readRDS("./gyn.combined_sub.rds")
Idents(gyn.combined_sub)="new_cell_type"
gyn.combined.core=subset(x = gyn.combined_sub, idents = c("mesophyll","BS","epidermis","guard","phloem","xylem","proliferating"))
DefaultAssay(gyn.combined.core)='RNA'
gyn.combined.core$new_cell_type=factor(gyn.combined.core$new_cell_type,
                                       levels=c('mesophyll','BS','xylem','phloem',
                                                'epidermis','guard','proliferating'))
colour_bar=(c("#238443", "#1d91c0",
              "#9970ab","#f1b6da",
              "#dfc27d", "#fec44f",
              "#fbb4ae"))
Idents(gyn.combined.core)="sample_type"
gyn.combined.core_nuc=subset(x = gyn.combined.core, idents = c("nuc"))
DimPlot(gyn.combined.core_nuc)
DefaultAssay(gyn.combined.core_nuc)='RNA'
Idents(gyn.combined.core_nuc)="new_cell_type"
DimPlot(gyn.combined.core_nuc,cols=colour_bar)
rm(gyn.combined_sub)
rm(gyn.combined.core)
#########
DefaultAssay(gyn.combined.core_nuc)='RNA'
CA="Gyn-10X-HIC-v5-00029450"
PEPC=c("Gyn-10X-HIC-v5-00018699")
PPDK=c("Gyn-10X-HIC-v5-00028147")
ALAAT1=c("Gyn-10X-HIC-v5-00022400")
BASS2=c("Gyn-10X-HIC-v5-00005585")
NAD_ME1=c('Gyn-10X-HIC-v5-00027505')
Rubsico=c("Gyn-10X-HIC-v5-00002736")
RCA=c('Gyn-10X-HIC-v5-00018250')
FBA1=c("Gyn-10X-HIC-v5-00021030")
TKL1=c("Gyn-10X-HIC-v5-00019024")
Core_C4=c(CA,PEPC,PPDK,ALAAT1,BASS2,NAD_ME1,Rubsico,RCA,FBA1,TKL1)
Stacked_VlnPlot(seurat_object = gyn.combined.core_nuc, features = Core_C4,colors_use=colour_bar,x_lab_rotate = TRUE)
