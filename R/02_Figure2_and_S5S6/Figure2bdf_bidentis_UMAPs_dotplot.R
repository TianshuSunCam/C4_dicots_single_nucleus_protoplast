library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
###bid
setwd("youtpayh")
bid.combined=readRDS("./bidentis_combined_integrated_05.rds")
###cluster annotation and dimplot
Idents(bid.combined)="integrated_cell_type"
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da",
             "#dfc27d", "#fec44f",
             "#fbb4ae","#bdbdbd")
DimPlot(bid.combined, reduction = "umap",label=F,group.by = "integrated_cell_type",cols=colour_bar)
###
#######cluster annotation
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
bid.combined.core=subset(x = bid.combined, idents = c("mesophyll","BS","epidermis","guard","phloem","xylem","proliferating"))
bid.combined.core=readRDS("./bidentis_combined_core.rds")
DefaultAssay(bid.combined.core)='RNA'
global_range <- c(-0.5,2.5)
bid_Dot=DotPlot(object =  bid.combined.core,
        features = rev(query_gene),scale.by = 'size',group.by = 'integrated_cell_type',  col.min = -0.5,
        col.max = 2.5,scale.min = 0, scale.max = 100)+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14,angle = 30,vjust=1,hjust = 1))& scale_colour_gradientn(colours = custom_magma, limits = global_range)
bid_Dot
#####read cell type marker genes and Arabidopsis orthologs
DefaultAssay(bidentis.combined)='RNA'
Idents(bidentis.combined) <- "integrated_cell_type"
bid_all_markers <- FindAllMarkers(object = bidentis.combined)
write.table(all.markers, 'bid_all_markers_cell_type.tab',col.names = T,row.names = T,quote = F)

bid_all_markers=read.table("./bid_all_markers_cell_type.tab")
bid_all_markers=bid_all_markers[bid_all_markers$avg_log2FC > 0, ]
ara_bid=read.table("./new_names_besthit_bidentis_arabidopsis.out",header = F)

bid_all_markers$ATG=ara_bid$V2[match(bid_all_markers$gene, ara_bid$V1)]
###library(DOSE)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(tidyverse)
library("org.At.tair.db")
meso_genes=bid_all_markers[which(bid_all_markers$cluster=='mesophyll'),]$ATG[1:200]
meso_GO = enrichGO(meso_genes, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,readable = TRUE)
head(meso_GO)
dotplot(meso_GO, showCategory=30) + ggtitle("dotplot for mesophyll")

####
BS_genes=bid_all_markers[which(bid_all_markers$cluster=='BS'),]$ATG[1:200]
BS_GO = enrichGO(BS_genes, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,readable = TRUE)
head(BS_GO)
dotplot(BS_GO, showCategory=30) + ggtitle("dotplot for BS")
###
xylem_genes=bid_all_markers[which(bid_all_markers$cluster=='xylem'),]$ATG[1:200]
xylem_GO = enrichGO(xylem_genes, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,readable = TRUE)
head(xylem_GO)
dotplot(xylem_GO, showCategory=30) + ggtitle("dotplot for xylem")###

phloem_genes=bid_all_markers[which(bid_all_markers$cluster=='phloem'),]$ATG[1:200]
epidermis_genes=bid_all_markers[which(bid_all_markers$cluster=='epidermis'),]$ATG[1:200]
guard_genes=bid_all_markers[which(bid_all_markers$cluster=='guard'),]$ATG[1:200]
proliferating_genes=bid_all_markers[which(bid_all_markers$cluster=='proliferating'),]$ATG[1:200]
###
query=bid_all_markers[which(bid_all_markers$cluster=='proliferating'),]$ATG[1:200]
query_GO = enrichGO(query, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,readable = TRUE)
head(query_GO)
dotplot(query_GO, showCategory=30) + ggtitle("dotplot for query")
###
enrich_gene_list_for_test <- list('mesophyll'=meso_genes,'BS'=BS_genes,
                                  'xylem'=xylem_genes,'phloem'=phloem_genes,
                                  'epidermis'=epidermis_genes,'guard'=guard_genes,
                                  'proliferating'=proliferating_genes)
bid_test_enrich <- compareCluster(geneCluster = enrich_gene_list_for_test,ont = "BP", 
                              OrgDb = org.At.tair.db, pAdjustMethod = "BH",keyType="TAIR",
                              pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                              fun =  enrichGO)
dotplot(bid_test_enrich, showCategory=30) + ggtitle("dotplot")+coord_flip()
bid_GO_for_plot=bid_test_enrich@compareClusterResult
Descripton_to_keep=c("photosynthesis","photosynthesis, light reaction","cellular response to oxygen levels",
                     "electron transport chain",
                     "response to wounding","response to acid chemical",
                     "photosynthesis, dark reaction","NAD(P)H dehydrogenase complex assembly",
                     "glycine metabolic process","serine family amino acid catabolic process",
                     "response to auxin","xylem development","inorganic anion transport",
                     "phloem development","stem cell population maintenance",
                     "regionalization","pattern specification process",
                     "cuticle development","wax biosynthetic process","fatty acid biosynthetic process",
                     "monocarboxylic acid biosynthetic process",
                     "response to carbon dioxide","stomatal complex development",
                     "cell cycle process","microtubule-based process")

head(bid_GO_for_plot)
bid_GO_for_plot=bid_test_enrich@compareClusterResult
bid_GO_for_plot=bid_GO_for_plot[which(bid_GO_for_plot$Description %in% Descripton_to_keep),]
bid_GO_for_plot$Description=factor(bid_GO_for_plot$Description,levels = rev(Descripton_to_keep))

bid_GO_for_plot$GeneRatio_numeric <- sapply(bid_GO_for_plot$GeneRatio, function(x) {
  ratio <- as.numeric(unlist(strsplit(x, "/")))  # Split and convert to numeric
  ratio[1] / ratio[2]  # Divide the two numbers
})

bid_GO=ggplot(bid_GO_for_plot, aes(x=Cluster, y=Description,color= p.adjust)) +
  geom_point(aes(size=GeneRatio_numeric))+scale_size(limits=c(0,0.5),breaks=c(0.1,0.2,0.3,0.4,0.5))+
  #geom_segment(aes(x = 0, xend = Cluster, y =Description, yend = Description),size = 1,color='grey')+
  scale_colour_gradientn(colours = c( "#440154FF" ,"#414487FF", "#2A788EFF" ,"#22A884FF" ,"#7AD151FF", "#FDE725FF"), limits = c(0,0.0056))+
  theme(text = element_text(size=16)) +
  theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1, vjust =1, face = "plain"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
gyn_GO+bid_GO
write.csv(bid_test_enrich@compareClusterResult,"bid_GO_compareClusterResult.csv")
