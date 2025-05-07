library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(ggpubr)

custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
bid.combined=readRDS("bidentis_combined_integrated_05.rds")

Idents(bid.combined) <- "sample_type"
DefaultAssay(bid.combined)="RNA"
options(future.globals.maxSize = 1000 * 1024^2)
#all.varied.re.markers <- FindMarkers(bid.combined, ident.1 = "proto", ident.2 = "nuc", verbose = T)
#write.table(all.varied.re.markers, 'bidentis_integrated_nuc_vs_proto_marker.tab',col.names = T,row.names = T,quote = F)

###bid
bid_nuc_proto_marker=read.table("/Users/tianshusun/Documents/work_2023/2023_method_paper/R/Code_for_figures/Figure3/bidentis_integrated_nuc_vs_proto_marker.tab")
###bid_nuc_proto_marker=read.table("/Users/tianshusun/Documents/work_2023/2023_method_paper/R/Code_for_figures/Figure3/all_markers_re_vaired.tab")
bid_nuc_proto_marker$gene=rownames(bid_nuc_proto_marker)
ara_bid=read.table("/Users/tianshusun/Documents/work_2023/2023_method_paper/R/Code_for_figures/Figure4/new_names_besthit_bidentis_arabidopsis.out",header = F)
bid_nuc_proto_marker$ATG=ara_bid$V2[match(bid_nuc_proto_marker$gene, ara_bid$V1)]
###marker gene
###library(DOSE)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(tidyverse)
library("org.At.tair.db")
proto_genes=bid_nuc_proto_marker[which(bid_nuc_proto_marker$avg_log2FC>0),]$ATG
proto_genes_GO = enrichGO(proto_genes, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,readable = TRUE)
head(proto_genes_GO)
dotplot(proto_genes_GO, showCategory=30) + ggtitle("dotplot for protoplast")

nuc_genes=bid_nuc_proto_marker[which(bid_nuc_proto_marker$avg_log2FC<0),]$ATG
nuc_genes_GO = enrichGO(nuc_genes, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,readable = TRUE)
head(nuc_genes_GO)
dotplot(nuc_genes_GO, showCategory=30) + ggtitle("dotplot for nuclei")
###
enrich_gene_list_for_test <- list('protoplast'=proto_genes,'nuclei'=nuc_genes)
bid_test_enrich <- compareCluster(geneCluster = enrich_gene_list_for_test,ont = "BP", 
                                  OrgDb = org.At.tair.db, pAdjustMethod = "BH",keyType="TAIR",
                                  pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                                  fun =  enrichGO)
dotplot(bid_test_enrich, showCategory=30) + ggtitle("dotplot")+coord_flip()
bid_GO_for_plot=bid_test_enrich@compareClusterResult
#####
Descripton_to_keep_proto=c("response to hypoxia","response to wounding","response to fungus",
                           "response to acid chemical",
                           "response to water deprivation","response to cold",
                           "response to jasmonic acid",
                           "cell death",
                           "response to salicylic acid","response to toxic substance")
Descripton_to_keep_nuc=c("gene silencing",
                         "ncRNA metabolic process","DNA repair","cell cycle process",
                         "microtubule-based process","mRNA metabolic process","chromatin organization","plastid organization",
                         "histone modification","nuclear export")
Descripton_to_keep=c(Descripton_to_keep_nuc,rev(Descripton_to_keep_proto))

bid_GO_for_plot=bid_test_enrich@compareClusterResult
bid_GO_for_plot=bid_GO_for_plot[which(bid_GO_for_plot$Description %in% Descripton_to_keep),]
bid_GO_for_plot$Description=factor(bid_GO_for_plot$Description,levels = c(Descripton_to_keep))

bid_GO_for_plot$GeneRatio_numeric <- sapply(bid_GO_for_plot$GeneRatio, function(x) {
  ratio <- as.numeric(unlist(strsplit(x, "/")))  # Split and convert to numeric
  ratio[1] / ratio[2]  # Divide the two numbers
})
ggplot(bid_GO_for_plot, aes(x=Description, y=Cluster,color= p.adjust)) +
  geom_point(aes(size=GeneRatio_numeric))+
#scale_size(limits=c(0,0.06),breaks=c(0.01,0.02,0.03,0.04,0.05,0.06))+
  #geom_segment(aes(x = Description, xend = Description, y = 0, yend = Cluster),size = 1,color='#d9d9d9')+
  scale_colour_gradientn(colours = c( "#440154FF" ,"#414487FF", "#2A788EFF" ,"#22A884FF" ,"#7AD151FF", "#FDE725FF"))+
  theme(text = element_text(size=16)) +theme(legend.position="top")+
  theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1, vjust =1, face = "plain"),
        panel.grid.major = element_line(color = "#969696", size = 0.5,linetype = 2), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

filtered_rows_proto <- bid_GO_for_plot[bid_GO_for_plot$Description %in% Descripton_to_keep_proto, "geneID"]
gene_ids_proto <- unlist(strsplit(filtered_rows_proto, "/"))
unique_gene_ids_proto <- unique(gene_ids_proto)
unique_gene_ids_proto_bid=bid_nuc_proto_marker[which(bid_nuc_proto_marker$ATG %in% unique_gene_ids_proto),]$gene
#write.csv(unique_gene_ids_proto_bid,'unique_gene_ids_proto_bid.txt')
unique_gene_ids_proto_bid=read.csv('unique_gene_ids_proto_bid.txt')
unique_gene_ids_proto_bid=unique_gene_ids_proto_bid$x
bid.combined_add = AddModuleScore(object = bid.combined,features = list(unique_gene_ids_proto_bid),name = 'stress')
DefaultAssay(bid.combined_add)="RNA"
Idents(bid.combined_add) <- "integrated_cell_type"
FeaturePlot(bid.combined_add,features = 'stress1',split.by ="sample_type",label=T)& scale_colour_gradientn(colours = custom_magma)& theme(legend.position = "right")
FeaturePlot(bid.combined_add,features = 'stress1',split.by ="sample_type",label=F, min.cutoff = -0.1, max.cutoff = 0.6)& 
  scale_colour_gradientn(colours = custom_magma,limits = c(-0.1,0.6))& theme(legend.position = "right")
Idents(bid.combined_add) <- "sample_type"
VlnPlot(bid.combined_add,features = 'stress1',pt.size = 0,cols=c("#737373","#addd8e"))
y_max <- max(FetchData(bid.combined_add, vars = 'stress1',), na.rm = TRUE)
VlnPlot(
  bid.combined_add,
  features = 'stress1',
  group.by = "sample_type",
  cols = c("#737373", "#addd8e"),
  pt.size = 0,
  y.max = y_max + 0.1 # Add buffer to avoid cutting off p-values
) + 
  stat_compare_means(comparisons = list(c("nuc", "proto")), label = "p.signif")


colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da", "#decbe4", "#de77ae",
             "#dfc27d", "#fec44f","#ffed6f", "#ffed6f",
             "#fbb4ae","#bdbdbd")
Idents(bid.combined_add) <- "sample_type"
VlnPlot(bid.combined_add,features = 'stress1',pt.size = 0,cols=c("#737373","#addd8e"))
Idents(bid.combined_add) <- "integrated_cell_type"
VlnPlot(bid.combined_add,features = 'stress1',pt.size = 0,cols=colour_bar)
VlnPlot(object = bid.combined_add, features = "stress1",  split.by = "sample_type", assay="integrated", pt.size=0,cols=c("#737373","#addd8e"))
library("Nebulosa")
plot_density(bid.combined_add, "stress1")



filtered_rows_nuc <- bid_GO_for_plot[bid_GO_for_plot$Description %in% Descripton_to_keep_nuc, "geneID"]
gene_ids_nuc <- unlist(strsplit(filtered_rows_nuc, "/"))
unique_gene_ids_nuc <- unique(gene_ids_nuc)
unique_gene_ids_nuc_bid=bid_nuc_proto_marker[which(bid_nuc_proto_marker$ATG %in% unique_gene_ids_nuc),]$gene
#write.csv(unique_gene_ids_nuc_bid,'unique_gene_ids_nuc_bid.txt')
unique_gene_ids_nuc_bid=read.csv('unique_gene_ids_nuc_bid.txt')
unique_gene_ids_nuc_bid=unique_gene_ids_nuc_bid$x

bid.combined_add = AddModuleScore(object = bid.combined,features = list(unique_gene_ids_nuc_bid),name = 'nuc')
DefaultAssay(bid.combined_add)="RNA"
Idents(bid.combined_add) <- "integrated_cell_type"
FeaturePlot(bid.combined_add,features = 'nuc1',split.by ="sample_type",label=T)& scale_colour_gradientn(colours = custom_magma)& theme(legend.position = "right")
FeaturePlot(bid.combined_add,features = 'nuc1',split.by ="sample_type",label=F, min.cutoff =-0.1, max.cutoff = 0.3)& scale_colour_gradientn(colours = custom_magma,limits = c(-0.1,0.3))& theme(legend.position = "right")
Idents(bid.combined_add) <- "sample_type"
y_max <- max(FetchData(bid.combined_add, vars = 'nuc1',), na.rm = TRUE)
VlnPlot(
  bid.combined_add,
  features = 'nuc1',
  group.by = "sample_type",
  cols = c("#737373", "#addd8e"),
  pt.size = 0,
  y.max = y_max + 0.1 # Add buffer to avoid cutting off p-values
) + 
  stat_compare_means(comparisons = list(c("nuc", "proto")), label = "p.signif")




Idents(bid.combined_add) <- "integrated_cell_type"
Idents(bid.combined_add) <- "sample_type"
VlnPlot(bid.combined_add,features = 'nuc1',pt.size = 0,cols=c("#737373","#addd8e"))
VlnPlot(object = bid.combined_add, features = "nuc1",  split.by = "sample_type", assay="integrated", pt.size=0,cols=c("#737373","#addd8e"))

colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da", "#decbe4", "#de77ae",
             "#dfc27d", "#fec44f","#ffed6f", "#ffed6f",
             "#fbb4ae","#bdbdbd")
Idents(bid.combined_add) <- "sample_type"
bid.combined_add_nuc=subset(bid.combined_add,idents = "nuc")
bid.combined_add_proto=subset(bid.combined_add,idents = "proto")
Idents(bid.combined_add_nuc) <- "integrated_cell_type"
a=VlnPlot(bid.combined_add_nuc,features = 'nuc1',pt.size = 0,cols=colour_bar)
Idents(bid.combined_add_proto) <- "integrated_cell_type"
b=VlnPlot(bid.combined_add_proto,features = 'nuc1',pt.size = 0,cols=colour_bar)



