library(dplyr)
library(ggplot2)
library(ggsankey)
library(Seurat)
library(patchwork)
library(viridis)

######
Orthotable=read.table("./Orthogroups_singlerow_newnames_16species.table")
Orthotable$V1=gsub(':','',Orthotable$V1)

bidentis.combined.core.nuc=readRDS("./bid.combined.core_nuc.rds")
DefaultAssay(bidentis.combined.core.nuc)='RNA'
bidentis.combined.core.nuc=subset(x = bidentis.combined.core.nuc, idents = c("mesophyll","BS","epidermis","guard","phloem","xylem"))
celltype_markers_bid_nuc= FindAllMarkers(object = bidentis.combined.core.nuc,only.pos = T)

celltype_markers_bid_nuc$ortho <- Orthotable$V1[match(celltype_markers_bid_nuc$gene, Orthotable$V2)]

write.csv(celltype_markers_bid_nuc,'celltype_markers_nuc_bid_combined_core.csv',row.names = T )
celltype_markers_bid_nuc2=read.csv('celltype_markers_nuc_bid_combined_core.csv',row.names = 1,header = T)
#celltype_markers_bid_nuc <- celltype_markers_bid_nuc[celltype_markers_bid_nuc$cluster != "proliferating", ]

celltype_markers_ara=read.csv('celltype_markers_ara_proto_core.csv',row.names = 1,header = T)

#######
####
celltype_markers_bid_nuc_selected <- celltype_markers_bid_nuc %>%
  group_by(cluster) %>%
  slice_head(n = 200) %>%
  ungroup()  # Optional: ungroup if you don't need the grouping afterwards

# View the first few rows of the selected dataframe
head(celltype_markers_bid_nuc_selected)
celltype_markers_bid_nuc_selected
celltype_markers_bid_nuc_selected[celltype_markers_bid_nuc_selected$cluster == "xylem", ]
####
celltype_markers_ara_selected <- celltype_markers_ara %>%
  group_by(cluster) %>%
  slice_head(n = 200) %>%
  ungroup()  # Optional: ungroup if you don't need the grouping afterwards

# View the first few rows of the selected dataframe
head(celltype_markers_ara_selected)
celltype_markers_ara_selected[celltype_markers_ara_selected$cluster == "mesophyll", ]
intersect(celltype_markers_bid_nuc_selected[celltype_markers_bid_nuc_selected$cluster == "xylem", ]$ortho,celltype_markers_ara_selected[celltype_markers_ara_selected$cluster == "mesophyll", ]$ortho)
# 
bid <- data.frame(cluster = celltype_markers_bid_nuc_selected$cluster, ortho = celltype_markers_bid_nuc_selected$ortho)
ara <- data.frame(cluster = celltype_markers_ara_selected$cluster, ortho = celltype_markers_ara_selected$ortho)
intersect_g_b=intersect(celltype_markers_bid_nuc_selected$ortho,celltype_markers_ara_selected$ortho)
bid=bid[bid$ortho%in% intersect_g_b, ]
ara=ara[ara$ortho%in% intersect_g_b, ]
# Merge bid and ara on ortho
merged_df <- merge(ara, bid, by = "ortho", all = TRUE, suffixes = c("_ara", "_bid"))
table(merged_df$cluster_ara,merged_df$cluster_bid)
table(merged_df$cluster_ara)
table(merged_df$cluster_bid)
# Create a long format dataframe suitable for ggsankey
df <- merged_df %>%
  make_long(cluster_ara,cluster_bid)

colour_bar=rev(c("#238443", "#1d91c0",
                 "#9970ab","#f1b6da",
                 "#dfc27d", "#fec44f"))

df=as.data.frame(df)


df$node=factor(df$node,levels=rev(c('mesophyll','BS','xylem','phloem',
                                    'epidermis','guard')))

df$next_node=factor(df$next_node,levels=rev(c('mesophyll','BS','xylem','phloem',
                                              'epidermis','guard')))

ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, 
               fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "white") +
  geom_sankey_label(size = 5, color = "white", fill = "gray50") +
  scale_fill_manual(values =colour_bar)+
  #scale_fill_viridis_d() +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",       
        plot.title = element_text(hjust = .5)) +
  #  ggtitle("unique cell type markers in Fb top 20% in Gy")
  ggtitle("arabidopsis to bidentis")

######
merged_df[which(merged_df$cluster_ara=='mesophyll'&merged_df$cluster_bid=='BS'),]

ortho_meso_to_BS=(merged_df[which(merged_df$cluster_ara=='mesophyll'&merged_df$cluster_bid=='BS'),])$ortho

ara_meso_to_BS=celltype_markers_ara[which(celltype_markers_ara$ortho %in% ortho_meso_to_BS),]
#######


