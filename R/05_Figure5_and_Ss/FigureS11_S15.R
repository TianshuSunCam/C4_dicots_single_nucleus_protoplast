library(ggtree)
library(dplyr)
library(scCustomize)
library(Seurat)
# Define a named list of tree files
tree_files <- list(
  SpeciesTree = "./SpeciesTree_rooted.txt",
  PPa1_5 = "./OG0000467_tree.txt",
  PPa6 = "./OG0007017_tree.txt",
  PEPC = "./OG0000724_tree.txt",
  AMK = "./OG0006927_tree.txt",
  CA = "./OG0000462_tree.txt",
  NADME1 = "./OG0010813_tree.txt",
  NADME2 = "./OG0006182_tree.txt",
  NADPMEs = "./OG0000536_tree.txt",
  PEPCK = "./OG0002656_tree.txt"
)

# Function to read and plot tree
plot_tree <- function(file, title = NULL) {
  tree <- read.tree(file)
  p <- ggtree(tree, ladderize = TRUE) + 
    geom_tiplab(size = 3, color = "black")
  if (!is.null(title)) p <- p + ggtitle(title)
  print(p)
}

# Loop through the list and plot each tree
for (name in names(tree_files)) {
  plot_tree(tree_files[[name]], title = name)
}

library(dplyr)
library(scCustomize)
#####
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da",
             "#dfc27d", "#fec44f",
             "#fbb4ae","#f2f2f2")
###
# Load Seurat objects
bid_obj <- readRDS("./bid.combined.core_nuc.rds")
gyn_obj <- readRDS('./gyn.combined.core_nuc.rds')

# Define color scales
global_range <- c(-1, 2.5)

# Custom plotting function for violin plots
plot_violin <- function(obj, features, title = NULL) {
  Stacked_VlnPlot(
    seurat_object = obj,
    features = features,
    colors_use = colour_bar,
    x_lab_rotate = TRUE,
    pt.size = 0.1
  ) + ggtitle(title)
}

#--------------------#
# DotPlot: PPa1/5 and PPa6
#--------------------#
Gyn_ppa1_5 <- c("Gyn-10X-HIC-v5-00017826","Gyn-10X-HIC-v5-00026500","Gyn-10X-HIC-v5-00028442","Gyn-10X-HIC-v5-00008684")
bid_ppa1_5 <- c("FBID-00041790","FBID-00041798","F-bidentis-27-000289","FBID-00014731","FBID-00045483","F-bidentis-S0-001593","FBID-00030810")
Gyn_ppa6   <- "Gyn-10X-HIC-v5-00016050"
bid_ppa6   <- c("F-bidentis-S0-001902","FBID-00001691")

dot_a <- DotPlot(gyn_obj, features = rev(c(Gyn_ppa1_5, Gyn_ppa6)), scale.min = 0, scale.max = 50) +
  coord_flip() & scale_colour_gradientn(colours = custom_magma, limits = global_range)

dot_b <- DotPlot(bid_obj, features = rev(c(bid_ppa1_5, bid_ppa6)), scale.min = 0, scale.max = 50) +
  coord_flip() & scale_colour_gradientn(colours = custom_magma, limits = global_range)

dot_a / dot_b  # 

#--------------------#
# Violin Plots for Gene Sets
#--------------------#

gene_sets <- list(
  PEPC = list(
    gyn = c("Gyn-10X-HIC-v5-00004255","Gyn-10X-HIC-v5-00007944","Gyn-10X-HIC-v5-00018699","Gyn-10X-HIC-v5-00018698","Gyn-10X-HIC-v5-00026241"),
    bid = c("FBID-00060088","FBID-00018755","FBID-00018760","FBID-00030156","FBID-00001952","F-bidentis-31-000425","F-bidentis-1-002572")
  ),
  AMK = list(
    gyn = c("Gyn-10X-HIC-v5-00002486","Gyn-10X-HIC-v5-00027122","Gyn-10X-HIC-v5-00027123","Gyn-10X-HIC-v5-00029799"),
    bid = c("FBID-00004673")
  ),
  CA = list(
    gyn = c("Gyn-10X-HIC-v5-00002957","Gyn-10X-HIC-v5-00029450","Gyn-10X-HIC-v5-00009823","Gyn-10X-HIC-v5-00022821"),
    bid = c("FBID-00057958","FBID-00012542","FBID-00026272","F-bidentis-22-001104","F-bidentis-22-001105","F-bidentis-22-000213","FBID-00033040","FBID-00057959","F-bidentis-21-000208")
  ),
  NADME1 = list(
    gyn = c("Gyn-10X-HIC-v5-00027505"),
    bid = c("FBID-00044033","FBID-00017393")
  ),
  NADME2 = list(
    gyn = c("Gyn-10X-HIC-v5-00003051"),
    bid = c("FBID-00006488","FBID-00006487","F-bidentis-30-000711")
  ),
  NADPME = list(
    gyn = c("Gyn-10X-HIC-v5-00016798","Gyn-10X-HIC-v5-00004493","Gyn-10X-HIC-v5-00028561"),
    bid = c("F-bidentis-27-000643","FBID-00013633","F-bidentis-14-000158","FBID-00003378","F-bidentis-30-000031","F-bidentis-19-000475")
  ),
  PEPCK = list(
    gyn = c("Gyn-10X-HIC-v5-00015955","Gyn-10X-HIC-v5-00020869"),
    bid = c("FBID-00019122")
  )
)


for (gene in names(gene_sets)) {
  print(plot_violin(gyn_obj, gene_sets[[gene]]$gyn, paste0("Gyn - ", gene)))
  print(plot_violin(bid_obj, gene_sets[[gene]]$bid, paste0("Bid - ", gene)))
}
