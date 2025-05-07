# Load the libraries
library(dplyr)
library(ggplot2)
library(ggsankey)
###read othotable

celltype_markers_gyn_nuc=read.csv('celltype_markers_nuc_gyn_combined_core.csv',row.names = 1,header = T)
celltype_markers_bid_nuc=read.csv('celltype_markers_nuc_bid_combined_core.csv',row.names = 1,header = T)

####
celltype_markers_gyn_nuc_selected <- celltype_markers_gyn_nuc %>%
  group_by(cluster) %>%
  slice_head(n = 200) %>%
  ungroup()  # Optional: ungroup if you don't need the grouping afterwards

# View the first few rows of the selected dataframe
head(celltype_markers_gyn_nuc_selected)
######slice_head(prop = 0.2) 
celltype_markers_bid_nuc_selected <- celltype_markers_bid_nuc %>%
  group_by(cluster) %>%
  slice_head(n = 200)%>%
  ungroup()  # Optional: ungroup if you don't need the grouping afterwards

# View the first few rows of the selected dataframe
head(celltype_markers_bid_nuc_selected)

# 
gyn <- data.frame(cluster = celltype_markers_gyn_nuc_selected$cluster, ortho = celltype_markers_gyn_nuc_selected$ortho)
bid <- data.frame(cluster = celltype_markers_bid_nuc_selected$cluster, ortho = celltype_markers_bid_nuc_selected$ortho)
intersect_g_b=intersect(celltype_markers_gyn_nuc_selected$ortho,celltype_markers_bid_nuc_selected$ortho)
gyn=gyn[gyn$ortho%in% intersect_g_b, ]
bid=bid[bid$ortho%in% intersect_g_b, ]

# Merge gyn and bid on ortho
merged_df <- merge(gyn, bid, by = "ortho", all = TRUE, suffixes = c("_gyn", "_bid"))
table(merged_df$cluster_gyn,merged_df$cluster_bid)
table(merged_df$cluster_gyn)
table(merged_df$cluster_bid)
# Create a long format dataframe suitable for ggsankey

df <- merged_df %>%
  make_long(cluster_gyn,cluster_bid)

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
  ggtitle("gynandra to bidentis")

####
library(dplyr)
gyn <- gyn %>%
  distinct(cluster, ortho, .keep_all = TRUE)
bid <- bid %>%
  distinct(cluster, ortho, .keep_all = TRUE)
merged_df <- merge(gyn, bid, by = "ortho", all = TRUE, suffixes = c("_gyn", "_bid"))

tair_desc=read.csv("Orthogroups_descriptions.tsv",header=F,sep = '\t',row.names = 1)

# Convert rownames of tair_desc into a column for joining
tair_desc <- tair_desc %>%
  tibble::rownames_to_column(var = "ortho")
# Merge the dataframes using left_join
merged_df <- merged_df %>%
  left_join(tair_desc, by = "ortho")
write.csv(merged_df,"4C_gyn_bid_ortho_description.csv")
