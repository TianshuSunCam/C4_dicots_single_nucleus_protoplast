library(UpSetR)

all.markers_bid <- read.table("./bid_all_markers_cell_type.tab")

all.markers_bid<-read.table("/Users/tianshusun/Documents/work_2023/2023_method_paper/R/Code_for_figures/Figure2/bid_all_markers_cell_type.tab")
all.markers_bid<-all.markers_bid[all.markers_bid$avg_log2FC > 0, ]
  
meso_celltype_marker=(all.markers_bid[which((all.markers_bid$cluster) =='mesophyll'), ])$gene[1:200]
BS_celltype_marker=(all.markers_bid[which((all.markers_bid$cluster) =='BS'), ]) $gene[1:200]
xylem_celltype_marker=(all.markers_bid[which((all.markers_bid$cluster) =='xylem'), ]) $gene[1:200]
phloem_celltype_marker=(all.markers_bid[which((all.markers_bid$cluster) =='phloem'), ]) $gene[1:200]
epidermis_celltype_marker=(all.markers_bid[which((all.markers_bid$cluster) =='epidermis'), ]) $gene[1:200]
guard_celltype_marker=(all.markers_bid[which((all.markers_bid$cluster) =='guard'), ]) $gene[1:200]
proliferating_celltype_marker=(all.markers_bid[which((all.markers_bid$cluster) =='proliferating'), ])$gene[1:200] 

listinput <- 
  list(mesophyll= c(as.character(meso_celltype_marker)), 
       BS = c(as.character(BS_celltype_marker)), 
       xylem = c(as.character(xylem_celltype_marker)),
       phloem= c(as.character(phloem_celltype_marker)), 
       epidermis= c(as.character(epidermis_celltype_marker)), 
       guard= c(as.character(guard_celltype_marker)),
       proliferating = c(as.character(proliferating_celltype_marker)))

set_order <- c("mesophyll", "BS", "xylem","phloem","epidermis","guard","proliferating")

upset(fromList(listinput),matrix.color = "#2171b5", main.bar.color = "#b3cde3",nsets = 6,
      sets = set_order,   keep.order = TRUE, 
      mainbar.y.label = "shared celltype marker genes", mainbar.y.max = NULL,
      sets.bar.color = "#b3cde3", sets.x.label = "No. of marker genes",
      text.scale = c(2, 2, 2, 2, 2, 2))
