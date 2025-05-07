library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(scCustomize)
library(ggpubr)
library(pheatmap)
library(stringi)
################read seurat object bidentis, nuclei only (subset from the integrated data)
bid.combined.core_nuc=readRDS("./bid.combined.core_nuc.rds")
Idents(bid.combined.core_nuc)="integrated_cell_type"
colour_bar=(c("#238443", "#1d91c0",
              "#9970ab","#f1b6da",
              "#dfc27d", "#fec44f",
              "#fbb4ae"))
DimPlot(bid.combined.core_nuc,cols=colour_bar)

####find mesophyll and BS markers
meso_BS.markers <- FindMarkers(bid.combined.core_nuc, ident.1 = 'mesophyll', ident.2 = 'BS')
meso_BS.markers_higher_meso=meso_BS.markers[which(meso_BS.markers$avg_log2FC>0),]
meso_BS.markers_higher_meso <- meso_BS.markers_higher_meso[order(meso_BS.markers_higher_meso$p_val_adj), ]
meso_BS.markers_higher_meso=meso_BS.markers_higher_meso[1:100,]
meso_BS.markers_higher_meso <- meso_BS.markers_higher_meso[order(meso_BS.markers_higher_meso$avg_log2FC), ]

meso_BS.markers_higher_BS=meso_BS.markers[which(meso_BS.markers$avg_log2FC<0),]
meso_BS.markers_higher_BS <- meso_BS.markers_higher_BS[order(meso_BS.markers_higher_BS$p_val_adj), ]
meso_BS.markers_higher_BS=meso_BS.markers_higher_BS[1:100,]
meso_BS.markers_higher_BS <- meso_BS.markers_higher_BS[order(meso_BS.markers_higher_BS$avg_log2FC), ]

qgene_meso=rownames(meso_BS.markers_higher_meso)
qgene_bs=rownames(meso_BS.markers_higher_BS)


Idents(bid.combined.core_nuc)='integrated_cell_type'

cluster.averages <- AverageExpression(bid.combined.core_nuc, return.seurat = T)
DoHeatmap(cluster.averages, features = qgene_bs, size = 3, draw.lines = FALSE)+  scale_fill_gradientn(colors = c("steelblue2", "white", "tomato2"))



a=pheatmap(cluster.averages[['RNA']][qgene_meso,],scale='row',cluster_cols = F,cluster_rows = T,
           border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
b=pheatmap(cluster.averages[['RNA']][qgene_bs,],scale='row',cluster_cols = F,cluster_rows = T,
           border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
order_row_meso =rev(a$tree_row$order) 
meso_data = data.frame(cluster.averages[['RNA']][qgene_meso,][order_row_meso,]) 

order_row_bs = b$tree_row$order
bs_data = data.frame(cluster.averages[['RNA']][qgene_bs,][order_row_bs,]) 

bid_data=rbind(meso_data,bs_data)

pheatmap(bid_data,scale='row',cluster_cols = F,cluster_rows = F,
         border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
####sort by the expression in mesophyll
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
meso_data = data.frame(cluster.averages[['RNA']][qgene_meso,])
meso_data[rev(order(meso_data$mesophyll)), ]
scale_data_meso=scale_rows(meso_data)
order_row_meso=rev(order(scale_data_meso$mesophyll))
scale_data_meso[order_row_meso, ]
pheatmap(scale_data_meso[order_row_meso, ],cluster_cols = F,cluster_rows = F,
          border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#####
bs_data = data.frame(cluster.averages[['RNA']][qgene_bs,])
bs_data[rev(order(bs_data$bsphyll)), ]
scale_data_bs=scale_rows(bs_data)
order_row_bs=rev(order(scale_data_bs$BS))
scale_data_bs[order_row_bs, ]
pheatmap(scale_data_bs[order_row_bs, ],cluster_cols = F,cluster_rows = F,
         border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
####
Breaks <- seq(min(min(scale_data_meso),min(scale_data_bs)),
              max(max(scale_data_meso),max(scale_data_bs)), length = 50)
Breaks <- seq(-2.5,2.5, length = 50)
pheatmap(scale_data_meso[order_row_meso, ],cluster_cols = F,cluster_rows = F,breaks = Breaks,
         border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(scale_data_bs[order_row_bs, ],cluster_cols = F,cluster_rows = F,breaks = Breaks,
         border_color=NA,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
####

