library(dplyr)
library(ggplot2)
library(ggsankey)
library(Seurat)
library(patchwork)
library(viridis)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

######
Orthotable=read.table("./Orthogroups_singlerow_newnames_16species.table")
Orthotable$V1=gsub(':','',Orthotable$V1)
gyn.combined_sub=readRDS("./gyn.combined_sub.rds")
Idents(gyn.combined_sub)="new_cell_type"
gyn.combined.core=subset(x = gyn.combined_sub, idents = c("mesophyll","BS","epidermis","guard","phloem","xylem"))
DefaultAssay(gyn.combined.core)='integrated'
gyn.combined.core$new_cell_type=factor(gyn.combined.core$new_cell_type, levels=c('mesophyll','BS','xylem','phloem', 'epidermis','guard'))
colour_bar=(c("#238443", "#1d91c0", "#9970ab","#f1b6da","#dfc27d", "#fec44f"))
DimPlot(gyn.combined.core,cols=colour_bar)
rm(gyn.combined_sub)
Idents(gyn.combined.core)="sample_type"
gyn.combined.core_nuc=subset(x = gyn.combined.core, idents = c("nuc"))
Idents(gyn.combined.core_nuc)="new_cell_type"
DefaultAssay(gyn.combined.core_nuc)='RNA'
celltype_markers_gyn_nuc= FindAllMarkers(object = gyn.combined.core_nuc,only.pos = T)
celltype_markers_gyn_nuc$ortho <- Orthotable$V1[match(celltype_markers_gyn_nuc$gene, Orthotable$V2)]
DimPlot(gyn.combined.core_nuc,cols=colour_bar)

celltype_markers_gyn_nuc$ortho <- Orthotable$V1[match(celltype_markers_gyn_nuc$gene, Orthotable$V2)]
write.csv(celltype_markers_gyn_nuc,'celltype_markers_nuc_gyn_combined_core.csv',row.names = T )
celltype_markers_gyn_nuc=read.csv('celltype_markers_nuc_gyn_combined_core.csv',row.names = 1,header = T)
celltype_markers_ara=read.csv('celltype_markers_ara_proto_core.csv',row.names = 1,header = T)

#######
####
celltype_markers_gyn_nuc_selected <- celltype_markers_gyn_nuc %>%
  group_by(cluster) %>%
  slice_head(n = 200) %>%
  ungroup()  # Optional: ungroup if you don't need the grouping afterwards

# View the first few rows of the selected dataframe
head(celltype_markers_gyn_nuc_selected)
####
celltype_markers_ara_selected <- celltype_markers_ara %>%
  group_by(cluster) %>%
  slice_head(n = 200) %>%
  ungroup()  # Optional: ungroup if you don't need the grouping afterwards

# View the first few rows of the selected dataframe
head(celltype_markers_ara_selected)

# 
gyn <- data.frame(cluster = celltype_markers_gyn_nuc_selected$cluster, ortho = celltype_markers_gyn_nuc_selected$ortho)
ara <- data.frame(cluster = celltype_markers_ara_selected$cluster, ortho = celltype_markers_ara_selected$ortho)
intersect_g_b=intersect(celltype_markers_gyn_nuc_selected$ortho,celltype_markers_ara_selected$ortho)
gyn=gyn[gyn$ortho%in% intersect_g_b, ]
ara=ara[ara$ortho%in% intersect_g_b, ]
# Merge gyn and ara on ortho
merged_df <- merge(ara, gyn, by = "ortho", all = TRUE, suffixes = c("_ara", "_gyn"))
table(merged_df$cluster_ara,merged_df$cluster_gyn)
table(merged_df$cluster_ara)
table(merged_df$cluster_gyn)
# Create a long format dataframe suitable for ggsankey
df <- merged_df %>%
  make_long(cluster_ara,cluster_gyn)

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
  ggtitle("arabidopsis to gynandra")

#####
merged_df[which(merged_df$cluster_ara=='mesophyll'&merged_df$cluster_gyn=='BS'),]

ortho_meso_to_BS=(merged_df[which(merged_df$cluster_ara=='mesophyll'&merged_df$cluster_gyn=='BS'),])$ortho

ara_meso_to_BS=celltype_markers_ara[which(celltype_markers_ara$ortho %in% ortho_meso_to_BS),]

#######
