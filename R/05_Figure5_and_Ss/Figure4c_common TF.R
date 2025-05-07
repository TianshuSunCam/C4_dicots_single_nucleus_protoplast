library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(scCustomize)
library(ggpubr)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

#####run the script for Figure 5a and b first
###intersect BS
BS_inter_ortho_list=intersect(gyn_marker_TF_bs$ortho ,bid_marker_TF_bs$ortho )

###intersect mesophyll
meso_inter_ortho_list=intersect(gyn_marker_TF_meso$ortho ,bid_marker_TF_meso$ortho )

###dotplot gynandra
qgene_BS_gyn=rownames(gyn_marker_TF_bs[which((gyn_marker_TF_bs$ortho) %in% BS_inter_ortho_list), ])
qgene_meso_gyn=rownames(gyn_marker_TF_meso[which((gyn_marker_TF_meso$ortho) %in% meso_inter_ortho_list), ])
gyn_marker_TF_bs[qgene_BS_gyn,]
table(gyn_marker_TF_bs[qgene_BS_gyn,]$ortho)
table(gyn_marker_TF_meso[qgene_meso_gyn,]$ortho)
##remove duplicated ortho and rearrange
gyn_meso_dataframe= gyn_marker_TF_meso[qgene_meso_gyn,]
gyn_meso_dataframe$gene=rownames(gyn_meso_dataframe)
# Load dplyr package
library(dplyr)
gyn_meso_dataframe <- gyn_meso_dataframe %>%
  group_by(ortho) %>%  # Group by 'ortho' column
  slice_min(p_val)  %>%     # Select the row with the minimum 'pvalue' in each group
  arrange(p_val)    # Select the row with the minimum 'pvalue' in each group
# View the result
head(gyn_meso_dataframe)
qgene_meso_gyn=gyn_meso_dataframe$gene

##remove duplicated ortho and rearrange
gyn_BS_dataframe= gyn_marker_TF_bs[qgene_BS_gyn,]
gyn_BS_dataframe$gene=rownames(gyn_BS_dataframe)
gyn_BS_dataframe <- gyn_BS_dataframe %>%
  group_by(ortho) %>%  # Group by 'ortho' column
  slice_min(p_val)  %>%     # Select the row with the minimum 'pvalue' in each group
  arrange(p_val)    # Select the row with the minimum 'pvalue' in each group
# View the result
head(gyn_BS_dataframe)
qgene_BS_gyn=gyn_BS_dataframe$gene
###########
qgene_gyn=c(qgene_BS_gyn,qgene_meso_gyn)
min(a$data$avg.exp.scaled)
max(a$data$avg.exp.scaled)
global_range <- c(-1.5,2.3)
a=DotPlot( gyn.combined.core_nuc, features =rev(qgene_gyn),scale.min = 0,
         scale.max = 50)+coord_flip() & scale_colour_gradientn(colours = custom_magma,limits = global_range)
gyn_marker_TF_bs[qgene_BS_gyn,]
gyn_marker_TF_meso[qgene_meso_gyn,]
######################
######################
######################
###dotplot bidentis
qgene_BS_bid= rownames(bid_marker_TF_bs[which((bid_marker_TF_bs$ortho) %in% BS_inter_ortho_list), ])
qgene_meso_bid= rownames(bid_marker_TF_meso[which((bid_marker_TF_meso$ortho) %in% meso_inter_ortho_list), ])
table(bid_marker_TF_bs[qgene_BS_bid,]$ortho)
table(bid_marker_TF_meso[qgene_meso_bid,]$ortho)
##remove duplicated ortho
bid_meso_dataframe= bid_marker_TF_meso[qgene_meso_bid,]
bid_meso_dataframe$gene=rownames(bid_meso_dataframe)
bid_meso_dataframe <- bid_meso_dataframe %>%
  group_by(ortho) %>%  # Group by 'ortho' column
  slice_min(p_val)  %>%     # Select the row with the minimum 'pvalue' in each group
  arrange(p_val)    # Select the row with the minimum 'pvalue' in each group

head(bid_meso_dataframe)
qgene_meso_bid=bid_meso_dataframe$gene
qgene_bid=c(qgene_BS_bid,qgene_meso_bid)
###order bidentis data as gynandra
order_data_gyn=rbind(gyn_marker_TF_meso[qgene_meso_gyn,],gyn_marker_TF_bs[qgene_BS_gyn,])
order_data_gyn$gene=rownames(order_data_gyn)
order_data_gyn$gene <- factor(order_data_gyn$gene, levels = qgene_gyn)
ordered_data_gyn <- order_data_gyn[order(order_data_gyn$gene), ]
ortho_order=ordered_data_gyn$ortho

order_data_bid=rbind(bid_marker_TF_meso[qgene_meso_bid,],bid_marker_TF_bs[qgene_BS_bid,])
order_data_bid$ortho <- factor(order_data_bid$ortho, levels = ortho_order)
order_data_bid <- order_data_bid[order(order_data_bid$ortho), ]
qgene_bid=rownames(order_data_bid)
min(b$data$avg.exp.scaled)
max(b$data$avg.exp.scaled)

b=DotPlot( bid.combined.core_nuc, features =rev(qgene_bid),  scale.min = 0,
         scale.max = 50)+coord_flip() & scale_colour_gradientn(colours = custom_magma,limits = global_range)

a+b
bid_marker_TF_bs[qgene_BS_bid,]
bid_marker_TF_meso[qgene_meso_bid,]