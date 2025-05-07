################load the required libraries################
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(scCustomize)
library(ggpubr)
library(pheatmap)
################read seurat object bidentis
bid.combined_sub=readRDS("./bidentis_combined_integrated_05.rds")
Idents(bid.combined_sub)="integrated_cell_type"
bid.combined.core=subset(x = bid.combined_sub, idents = c("mesophyll","BS","epidermis","guard","phloem","xylem"))
DefaultAssay(bid.combined.core)='integrated_cell_type'
bid.combined.core$integrated_cell_type=factor(bid.combined.core$integrated_cell_type, levels=c('mesophyll','BS','xylem','phloem', 'epidermis','guard'))
colour_bar=(c("#238443", "#1d91c0", "#9970ab","#f1b6da","#dfc27d", "#fec44f"))
DimPlot(bid.combined.core,cols=colour_bar)
rm(bid.combined_sub)
Idents(bid.combined.core)="sample_type"
bid.combined.core_nuc=subset(x = bid.combined.core, idents = c("nuc"))
Idents(bid.combined.core_nuc)="integrated_cell_type"
DefaultAssay(bid.combined.core_nuc)='RNA'

####identify marker genes between mesophyll and BS
meso_BS.markers_bid <- FindMarkers(bid.combined.core_nuc, ident.1 = 'mesophyll', ident.2 = 'BS')
meso_BS.markers_bid_hgiher_meso=meso_BS.markers_bid[which(meso_BS.markers_bid$avg_log2FC>0),]
meso_BS.markers_bid_hgiher_BS=meso_BS.markers_bid[which(meso_BS.markers_bid$avg_log2FC<0),]

####read the orthogruops for bidentis, gynandra and Arabidopsis
Orthotable=read.table("./Orthogroups_singlerow_newnames_16species.table")
Orthotable$V1=gsub(':',"",Orthotable$V1)
rownames(Orthotable)=Orthotable$V2
#######match marker with orthogruop
meso_BS.markers_bid$ortho=Orthotable[rownames(meso_BS.markers_bid),]$V1

#####read Arabidopsis thaliana Transcription Factors downloaded from https://planttfdb.gao-lab.org/
TF_type_table=read.table("Ath_TF_list.txt", row.names = NULL,header=T)
colnames(TF_type_table)=c('TF_ID','TF.Locus.ID','TF.Family.Name')
TF_type_table <- TF_type_table %>%
  filter(duplicated(TF.Locus.ID) == FALSE)
rownames(TF_type_table)=TF_type_table$TF.Locus.ID
################match TF family with orthogruop################
AtTF_Orthotable=Orthotable[which(Orthotable$V2 %in% TF_type_table$TF.Locus.ID),]
AtTF_Orthotable$type=TF_type_table[AtTF_Orthotable$V2,]$TF.Family.Name 
################find the best hit in arabidopsis for bidentis genes################
ara_bid=read.table("./new_names_besthit_bidentis_arabidopsis.out",header = F)
################find marker that is a TF################
meso_BS.markers_bid$ATG=ara_bid$V2[match(rownames(meso_BS.markers_bid), ara_bid$V1)]
bid_marker_TF=meso_BS.markers_bid[which(meso_BS.markers_bid$ATG %in% TF_type_table$TF.Locus.ID),]
################add TF family
bid_marker_TF$type <- AtTF_Orthotable$type[match(bid_marker_TF$ATG, AtTF_Orthotable$V2)]
freq_table <- table(bid_marker_TF$type)
bid_marker_TF$freq <- freq_table[bid_marker_TF$type]
bid_marker_TF_sorted <- bid_marker_TF[order(-bid_marker_TF$p_val_adj), ]
#bid_marker_TF_sorted <- bid_marker_TF_sorted[order(-bid_marker_TF_sorted$freq), ]
bid_marker_TF_sorted <- bid_marker_TF_sorted[order(bid_marker_TF_sorted$avg_log2FC), ]
qgene_bid=rownames(bid_marker_TF_sorted)
qgene_meso_bid=rownames(bid_marker_TF_sorted[which(bid_marker_TF_sorted$avg_log2FC>0),])
qgene_bs_bid=rownames(bid_marker_TF_sorted[which(bid_marker_TF_sorted$avg_log2FC<0),])
cluster.averages_bid <- AverageExpression(bid.combined.core_nuc, return.seurat = TRUE)
###
a_bid=pheatmap(cluster.averages_bid[['RNA']][qgene_meso_bid,],scale='row',cluster_cols = F,cluster_rows = T,
           border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
b_bid=pheatmap(cluster.averages_bid[['RNA']][qgene_bs_bid,],scale='row',cluster_cols = F,cluster_rows = T,
           border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
order_row_meso_bid = a_bid$tree_row$order
meso_TF_data_bid = data.frame(cluster.averages_bid[['RNA']][qgene_meso_bid,][order_row_meso_bid,]) 
order_row_bs_bid = b_bid$tree_row$order
bs_TF_data_bid = data.frame(cluster.averages_bid[['RNA']][qgene_bs_bid,][order_row_bs_bid,]) 
bid_TF_data=rbind(meso_TF_data_bid,bs_TF_data_bid)

####
####sort by the expression in mesophyll
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
meso_TF_data_bid = data.frame(cluster.averages_bid[['RNA']][qgene_meso_bid,])
meso_TF_data_bid[rev(order(meso_TF_data_bid$mesophyll)), ]
scale_data_meso_TF_bid=scale_rows(meso_TF_data_bid)
order_row_meso_TF_bid=rev(order(scale_data_meso_TF_bid$mesophyll))
scale_data_meso_TF_bid[order_row_meso_TF_bid, ]

#####
bs_TF_data_bid = data.frame(cluster.averages_bid[['RNA']][qgene_bs_bid,])
bs_TF_data_bid[rev(order(bs_TF_data_bid$bsphyll)), ]
scale_data_bs_TF_bid=scale_rows(bs_TF_data_bid)
order_row_bs_TF_bid=rev(order(scale_data_bs_TF_bid$BS))
scale_data_bs_TF_bid[order_row_bs_TF_bid, ]

####
Breaks <- seq(min(min(scale_data_meso_TF_bid),min(scale_data_bs_TF_bid)),
              max(max(scale_data_meso_TF_bid),max(scale_data_bs_TF_bid)), length = 50)
Breaks <- seq(-2,2.5, length = 50)
pheatmap(scale_data_meso_TF_bid[order_row_meso_TF_bid, ],cluster_cols = F,cluster_rows = F,breaks = Breaks,
         border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(scale_data_bs_TF_bid[order_row_bs_TF_bid, ],cluster_cols = F,cluster_rows = F,,breaks = Breaks,
         border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

###
bid_marker_TF_meso=bid_marker_TF_sorted[which(bid_marker_TF_sorted$avg_log2FC>0),]
bid_marker_TF_meso_type=as.data.frame(table(bid_marker_TF_meso$type))
bid_marker_TF_meso_type$color=rep('meso',nrow(bid_marker_TF_meso_type))

bid_marker_TF_bs=bid_marker_TF_sorted[which(bid_marker_TF_sorted$avg_log2FC<0),]
bid_marker_TF_bs_type=as.data.frame(table(bid_marker_TF_bs$type))
bid_marker_TF_bs_type$color=rep('bs',nrow(bid_marker_TF_bs_type))

bid_marker_TF_type=rbind(bid_marker_TF_bs_type,bid_marker_TF_meso_type)
colnames(bid_marker_TF_type)=c("TF",'Freq','color')

library(tidyr)
library(dplyr)
df_complete_bid <- bid_marker_TF_type %>%
  complete(TF, color, fill = list(Freq = 0))

# Assuming your dataframe is named 'df'

# Create separate columns for 'bs' and 'meso' Freq values
df_sorted_bid <- df_complete_bid %>%
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
df_sorted_bid <- df_sorted_bid %>%
  mutate(TF = factor(TF, levels =rev(unique(TF)) ))

df_sorted_bid$color=factor(df_sorted_bid$color,levels = c('meso','bs'))

