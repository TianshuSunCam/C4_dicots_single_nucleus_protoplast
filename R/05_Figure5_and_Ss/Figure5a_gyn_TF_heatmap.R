library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(scCustomize)
library(ggpubr)
library(pheatmap)
library(stringi)


custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

#######subset nuclei data of major cluster from integrate seurat obj
gyn.combined_sub=readRDS("gyn.combined_sub.rds")
Idents(gyn.combined_sub)="new_cell_type"
DefaultAssay(gyn.combined_sub)='RNA'
gyn.combined.core=subset(x = gyn.combined_sub, idents = c("mesophyll","BS","epidermis","guard","phloem","xylem","proliferating"))
Idents(gyn.combined.core)="sample_type"
gyn.combined.core_nuc=subset(x = gyn.combined.core, idents = c("nuc"))
Idents(gyn.combined.core_nuc)="new_cell_type"
DimPlot(gyn.combined.core_nuc,cols=colour_bar)
rm(gyn.combined_sub)
rm(gyn.combined.core)
DefaultAssay(gyn.combined.core_nuc)='RNA'
####find mesophyll and BS markers
meso_BS.markers <- FindMarkers(gyn.combined.core_nuc, ident.1 = 'mesophyll', ident.2 = 'BS')
meso_BS.markers_hgiher_meso=meso_BS.markers[which(meso_BS.markers$avg_log2FC>0),]
meso_BS.markers_hgiher_BS=meso_BS.markers[which(meso_BS.markers$avg_log2FC<0),]

#######reading orthogroup
Orthotable=read.table("Orthogroups_singlerow_newnames_16species.table")
Orthotable$V1=gsub(':',"",Orthotable$V1)
rownames(Orthotable)=Orthotable$V2
#######match marker with orthogruop
meso_BS.markers$ortho=Orthotable[rownames(meso_BS.markers),]$V1


######read Arabidopsis thaliana Transcription Factors downloaded from https://planttfdb.gao-lab.org/
TF_type_table=read.table("Ath_TF_list.txt", row.names = NULL,header=T)
colnames(TF_type_table)=c('TF_ID','TF.Locus.ID','TF.Family.Name')
TF_type_table <- TF_type_table %>%
  filter(duplicated(TF.Locus.ID) == FALSE)
rownames(TF_type_table)=TF_type_table$TF.Locus.ID

###match TF family with orthogruop
AtTF_Orthotable=Orthotable[which(Orthotable$V2 %in% TF_type_table$TF.Locus.ID),]
AtTF_Orthotable$type=TF_type_table[AtTF_Orthotable$V2,]$TF.Family.Name 
###find the best hit in arabidopsis for gynandra genes
ara_gyn=read.table("besthit_gyn_to_arabidopsis.out",header = F)
ara_gyn$V1=sapply(stri_split_fixed(ara_gyn$V1, '.', n = -1, simplify = FALSE), `[`, 1)
ara_gyn$V2=sapply(stri_split_fixed(ara_gyn$V2, '.', n = -1, simplify = FALSE), `[`, 1)
###find marker that is a TF
meso_BS.markers$ATG=ara_gyn$V2[match(rownames(meso_BS.markers), ara_gyn$V1)]
gyn_marker_TF=meso_BS.markers[which(meso_BS.markers$ATG %in% TF_type_table$TF.Locus.ID),]
#####add TF family
gyn_marker_TF$type <- AtTF_Orthotable$type[match(gyn_marker_TF$ATG, AtTF_Orthotable$V2)]

freq_table <- table(gyn_marker_TF$type)
gyn_marker_TF$freq <- freq_table[gyn_marker_TF$type]

gyn_marker_TF_sorted <- gyn_marker_TF[order(-gyn_marker_TF$p_val_adj), ]
#gyn_marker_TF_sorted <- gyn_marker_TF_sorted[order(-gyn_marker_TF_sorted$freq), ]
gyn_marker_TF_sorted <- gyn_marker_TF_sorted[order(gyn_marker_TF_sorted$avg_log2FC), ]

qgene=rownames(gyn_marker_TF_sorted)
qgene_meso=rownames(gyn_marker_TF_sorted[which(gyn_marker_TF_sorted$avg_log2FC>0),])
qgene_bs=rownames(gyn_marker_TF_sorted[which(gyn_marker_TF_sorted$avg_log2FC<0),])


Idents(gyn.combined.core_nuc)='new_cell_type'
cluster.averages <- AverageExpression(gyn.combined.core_nuc, return.seurat = T)



a=pheatmap(cluster.averages[['RNA']][qgene_meso,],scale='row',cluster_cols = F,cluster_rows = T,
           border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
b=pheatmap(cluster.averages[['RNA']][qgene_bs,],scale='row',cluster_cols = F,cluster_rows = T,
           border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
order_row_meso = a$tree_row$order
meso_TF_data = data.frame(cluster.averages[['RNA']][qgene_meso,][order_row_meso,]) 

order_row_bs = b$tree_row$order
bs_TF_data = data.frame(cluster.averages[['RNA']][qgene_bs,][order_row_bs,]) 

gyn_TF_data=rbind(meso_TF_data,bs_TF_data)


####sort by the expression in mesophyll
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
meso_TF_data = data.frame(cluster.averages[['RNA']][qgene_meso,])
meso_TF_data[rev(order(meso_TF_data$mesophyll)), ]
scale_data_meso_TF=scale_rows(meso_TF_data)
order_row_meso_TF=rev(order(scale_data_meso_TF$mesophyll))
scale_data_meso_TF[order_row_meso_TF, ]

#####
bs_TF_data = data.frame(cluster.averages[['RNA']][qgene_bs,])
bs_TF_data[rev(order(bs_TF_data$bsphyll)), ]
scale_data_bs_TF=scale_rows(bs_TF_data)
order_row_bs_TF=rev(order(scale_data_bs_TF$BS))
scale_data_bs_TF[order_row_bs_TF, ]

####
Breaks <- seq(min(min(scale_data_meso_TF),min(scale_data_bs_TF)),
              max(max(scale_data_meso_TF),max(scale_data_bs_TF)), length = 50)
Breaks <- seq(-2,2.5, length = 50)
pheatmap(scale_data_meso_TF[order_row_meso_TF, ],cluster_cols = F,cluster_rows = F,breaks = Breaks,
         border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(scale_data_bs_TF[order_row_bs_TF, ],cluster_cols = F,cluster_rows = F,breaks = Breaks,
         border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
####
gyn_marker_TF_meso=gyn_marker_TF_sorted[which(gyn_marker_TF_sorted$avg_log2FC>0),]
gyn_marker_TF_meso_type=as.data.frame(table(gyn_marker_TF_meso$type))
gyn_marker_TF_meso_type$color=rep('meso',nrow(gyn_marker_TF_meso_type))

gyn_marker_TF_bs=gyn_marker_TF_sorted[which(gyn_marker_TF_sorted$avg_log2FC<0),]
gyn_marker_TF_bs_type=as.data.frame(table(gyn_marker_TF_bs$type))
gyn_marker_TF_bs_type$color=rep('bs',nrow(gyn_marker_TF_bs_type))

gyn_marker_TF_type=rbind(gyn_marker_TF_bs_type,gyn_marker_TF_meso_type)
colnames(gyn_marker_TF_type)=c("TF",'Freq','color')

library(tidyr)
library(dplyr)
df_complete_gyn <- gyn_marker_TF_type %>%
  complete(TF, color, fill = list(Freq = 0))

# Assuming your dataframe is named 'df'

# Create separate columns for 'bs' and 'meso' Freq values
df_sorted_gyn <- df_complete_gyn %>%
  # Get the Freq where color is 'bs'
  group_by(TF) %>%
  mutate(Freq_bs = ifelse(color == 'bs', Freq, NA)) %>%
  fill(Freq_bs, .direction = 'updown') %>%
  # Get the Freq where color is 'meso'
  mutate(Freq_meso = ifelse(color == 'meso', Freq, NA)) %>%
  fill(Freq_meso, .direction = 'updown') %>%
  ungroup() %>%
  # Sort by Freq_bs (desc), then by Freq_meso (desc)
  arrange(desc(Freq_bs), desc(Freq_meso))
# Convert TF to a factor with levels in the current order
df_sorted_gyn <- df_sorted_gyn %>%
  mutate(TF = factor(TF, levels =rev(unique(TF)) ))

df_sorted_gyn$color=factor(df_sorted_gyn$color,levels = c('meso','bs'))
