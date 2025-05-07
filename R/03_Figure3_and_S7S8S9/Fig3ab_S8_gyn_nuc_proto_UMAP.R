library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
###gyn
gyn.combined=readRDS("./gyn.combined_sub.rds")
DefaultAssay(gyn.combined)='RNA'
###UMAP dimplot
colour_bar=c("#737373","#addd8e")
DimPlot(gyn.combined, reduction = "umap",group.by ='sample_type',label=F,cols=colour_bar)

###cluster annotation and dimplot
Idents(gyn.combined)="new_cell_type"
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da", "#decbe4", "#de77ae",
             "#dfc27d", "#fec44f","#ffed6f", "#ffed6f",
             "#fbb4ae","#bdbdbd")
DimPlot(gyn.combined, reduction = "umap",label=F,group.by = "new_cell_type",cols=colour_bar)
####
PPDK=c("Gyn-10X-HIC-v5-00028147")
PEPC=c("Gyn-10X-HIC-v5-00018699")
NAD_ME1=c('Gyn-10X-HIC-v5-00027505')
TKL1=c("Gyn-10X-HIC-v5-00019024")
ACL5=c('Gyn-10X-HIC-v5-00016900')
TMO6='Gyn-10X-HIC-v5-00017445'
KCS10=c('Gyn-10X-HIC-v5-00028867')
FAMA=c('Gyn-10X-HIC-v5-00027865')
HTA2=c('Gyn-10X-HIC-v5-00008303')
POK2='Gyn-10X-HIC-v5-00005859'

query_gene=c(PPDK,PEPC,NAD_ME1,TKL1,
             ACL5,TMO6,KCS10,FAMA,
             HTA2,POK2)
Idents(gyn.combined)="new_cell_type"
gyn.combined.core=subset(x = gyn.combined, idents = c("mesophyll","BS","epidermis","guard","phloem","xylem","proliferating"))
DefaultAssay(gyn.combined.core)='RNA'
DotPlot(gyn.combined.core, features = query_gene, split.by = "sample_type",cols = c("#000000","#41ab5d"))

####
gyn_celltype_composition=table(gyn.combined@meta.data$sample_type,gyn.combined@meta.data$new_cell_type)
gyn_celltype_composition=gyn_celltype_composition[,1:11]
percentages_gyn <- t(apply(gyn_celltype_composition, 1, function(x) (x / sum(x)) * 100))
# Reshape the matrix into a long format
library(reshape2)
percentages_gyn_long <- melt(percentages_gyn, varnames = c("Type", "Cell_Type"), 
                  value.name = "Percentage")
percentages_gyn_long$Cell_Type=factor(percentages_gyn_long$Cell_Type,
                                      levels=rev(c("mesophyll","BS","xylem","phloem",
                                                   "procambium","companion-SE",
                                                   "epidermis","guard",
                                                   "trichome-myb","trichome-sty",
                                                   "proliferating","ua")))
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da", "#decbe4", "#de77ae",
             "#dfc27d", "#fec44f","#ffed6f", "#ffed6f",
             "#fbb4ae","#bdbdbd")
ggplot(percentages_gyn_long,aes(x=Type,y=Percentage,fill=Cell_Type))+
  geom_bar(stat='identity',color='gray')+
  geom_text(aes(label=round(Percentage,0)), size=4, position=position_stack(vjust=.5))+
  scale_fill_manual(values=rev(colour_bar[1:11]))+
 #scale_y_continuous(labels = scales::percent_format())+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
#######

