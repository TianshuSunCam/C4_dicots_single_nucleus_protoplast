library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
###bid
bid.combined=readRDS("./bidentis_combined_integrated_05.rds")

###UMAP dimplot
colour_bar=c("#737373","#addd8e")
DimPlot(bid.combined, reduction = "umap",group.by ='sample_type',label=F,cols=colour_bar)

###
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da",
             "#dfc27d", "#fec44f",
             "#fbb4ae","#bdbdbd")
DimPlot(bid.combined, reduction = "umap",label=F,group.by = "integrated_cell_type",cols=colour_bar)
####
PPDK=c('F-bidentis-S0-002372')
PEPC='FBID-00060088' #'Fbid_CUFF.4218'
NADPME=c('FBID-00003378')
TKL1=c('FBID-00043968')
ACL5=c('FBID-00053265')
TMO6=c('F-bidentis-4-001300')
OPS=c('F-bidentis-S0-004020','F-bidentis-S0-001349')# Query: AT3G09070 'Gyn-10X-HIC-v5-00011689',
KCS10=c('FBID-00008631')
FAMA=c('FBID-00020874')
HTA2=c('F-bidentis-S0-000322')
POK2='FBID-00008020'#AT3G19050.1

query_gene=c(PPDK,PEPC,NADPME,TKL1,
             ACL5,TMO6,KCS10,FAMA,
             HTA2,POK2)
Idents(bid.combined)="integrated_cell_type"
DotPlot(bid.combined, features = query_gene, split.by = "sample_type",cols = c("#000000","#41ab5d"))
DotPlot(object =  bid.combined, features = c(query_gene),
        scale.by = 'size',split.by = 'sample_type',
        group.by = 'integrated_cell_type',cols = c("blue", "red"),
        scale.min = 0, scale.max = 100)+ 
                theme(axis.text.x = element_text(size = 14,angle = 30,vjust=1,hjust = 1))
####
bid_celltype_composition=table(bid.combined@meta.data$sample_type,bid.combined@meta.data$integrated_cell_type)
bid_celltype_composition=bid_celltype_composition[,1:7]
percentages_bid <- t(apply(bid_celltype_composition, 1, function(x) (x / sum(x)) * 100))
# Reshape the matrix into a long format
library(reshape2)
percentages_bid_long <- melt(percentages_bid, varnames = c("Type", "Cell_Type"), 
                  value.name = "Percentage")
percentages_bid_long$Cell_Type=factor(percentages_bid_long$Cell_Type,
                                      levels=rev(c("mesophyll","BS","xylem","phloem",
                                                   "epidermis","guard","proliferating","ua")))

ggplot(percentages_bid_long,aes(x=Type,y=Percentage,fill=Cell_Type))+
  geom_bar(stat='identity',color='gray')+
  geom_text(aes(label=round(Percentage,0)), size=4, position=position_stack(vjust=.5))+
  scale_fill_manual(values=rev(colour_bar[1:7]))+
 #scale_y_continuous(labels = scales::percent_format())+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
#######
Idents(bid.combined)='sample_type'
bid.combined.nuc=subset(x = bid.combined, idents = "nuc")
