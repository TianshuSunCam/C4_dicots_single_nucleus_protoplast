library(ggplot2)
library(Seurat)
#######
bidentis.combined=readRDS("./bidentis_combined_integrated_05.rds")
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da",
             "#dfc27d", "#fec44f",
             "#fbb4ae","#f2f2f2")
DimPlot(bidentis.combined, reduction = "umap",group.by ='integrated_cell_type',label=T,cols=colour_bar)
Idents(bidentis.combined)="integrated_cell_type"
bid.combined.core=subset(x = bidentis.combined, idents = c("mesophyll","BS","epidermis","guard","phloem","xylem","proliferating"))
bid.combined.core$new_cell_type=factor(bid.combined.core$new_cell_type,
                                       levels=c('mesophyll','BS','xylem','phloem',
                                                'epidermis','guard','proliferating'))
Idents(bid.combined.core)="sample_type"
bid.combined.core_nuc=subset(x = bid.combined.core, idents = c("nuc"))
DimPlot(bid.combined.core_nuc)
DefaultAssay(bid.combined.core_nuc)='RNA'
Idents(bid.combined.core_nuc)="integrated_cell_type"
DimPlot(bid.combined.core_nuc,cols=colour_bar)

rm(bid.combined_sub)
rm(bid.combined.core)
DefaultAssay(bid.combined.core_nuc)='RNA'
######################
bid.combined.core_nuc=readRDS("./bid.combined.core_nuc.rds")
Idents(bid.combined.core_nuc)="integrated_cell_type"
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da",
             "#dfc27d", "#fec44f",
             "#fbb4ae","#f2f2f2")

CA=c('F-bidentis-22-001105')
PEPC=c("FBID-00018755")#AT1G53310.3
PCK=c("FBID-00019122")
PPDK=c("F-bidentis-S0-002372")
RP1="F-bidentis-20-000141"
ALAAT1=c("FBID-00054102")
DiT1=c("FBID-00017792")
BASS2=c("FBID-00008921")
NADPME1=c('FBID-00003378')
mMDH1	=c('F-bidentis-44-000179')
Rubsico=c("FBID-00024335")
RCA=c('F-bidentis-51-000137')
Core_C4=c(CA,PEPC,PPDK,ALAAT1,BASS2,NADPME1,Rubsico,RCA,FBA1,TKL1)
Stacked_VlnPlot(seurat_object = bid.combined.core_nuc, features = Core_C4,colors_use=colour_bar,x_lab_rotate = F)
