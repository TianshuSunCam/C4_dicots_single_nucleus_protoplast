library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(ggpubr)

custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
gyn.combined=readRDS("./gyn.combined_sub.rds")

Idents(gyn.combined) <- "sample_type"
DefaultAssay(gyn.combined)="RNA"
options(future.globals.maxSize = 1000 * 1024^2)
#all.varied.re.markers <- FindMarkers(gyn.combined, ident.1 = "proto", ident.2 = "nuc", verbose = T)
#write.table(all.varied.re.markers, 'gynentis_integrated_nuc_vs_proto_marker.tab',col.names = T,row.names = T,quote = F)

###gyn
gyn_nuc_proto_marker=read.table(".gynandra_integrated_nuc_vs_proto_marker.tab")
###gyn_nuc_proto_marker=read.table("./all_markers_re_vaired.tab")
gyn_nuc_proto_marker$gene=rownames(gyn_nuc_proto_marker)
ara_gyn=read.table("./besthit_gyn_to_arabidopsis.out",header = F)
ara_gyn$V1=gsub("\\.mRNA[0-9]+","",ara_gyn$V1)
ara_gyn$V2=gsub("\\.[0-9]+","",ara_gyn$V2)
gyn_nuc_proto_marker$ATG=ara_gyn$V2[match(gyn_nuc_proto_marker$gene, ara_gyn$V1)]
###marker gene
###library(DOSE)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(tidyverse)
library("org.At.tair.db")
proto_genes=gyn_nuc_proto_marker[which(gyn_nuc_proto_marker$avg_log2FC>0),]$ATG
nuc_genes=gyn_nuc_proto_marker[which(gyn_nuc_proto_marker$avg_log2FC<0),]$ATG

###
enrich_gene_list_for_test <- list('protoplast'=proto_genes,'nuclei'=nuc_genes)
gyn_test_enrich <- compareCluster(geneCluster = enrich_gene_list_for_test,ont = "BP", 
                                  OrgDb = org.At.tair.db, pAdjustMethod = "BH",keyType="TAIR",
                                  pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                                  fun =  enrichGO)
dotplot(gyn_test_enrich, showCategory=30) + ggtitle("dotplot")+coord_flip()
gyn_GO_for_plot=gyn_test_enrich@compareClusterResult
#####
Descripton_to_keep_proto=c("response to hypoxia","response to water deprivation",
                           "response to acid chemical","response to wounding",
                           "response to fungus","response to salt stress", "toxin metabolic process",
                           "response to jasmonic acid","cell death",
                           "response to salicylic acid")
Descripton_to_keep_nuc=c("ncRNA metabolic process","plastid organization","mRNA metabolic process","gene silencing",
                         "chromatin organization",
                         "DNA repair","histone modification","cell cycle process",
                         "microtubule-based process", "nuclear export")

filtered_rows_proto <- gyn_GO_for_plot[gyn_GO_for_plot$Description %in% Descripton_to_keep_proto, "geneID"]
gene_ids_proto <- unlist(strsplit(filtered_rows_proto, "/"))
unique_gene_ids_proto <- unique(gene_ids_proto)
unique_gene_ids_proto_gyn=gyn_nuc_proto_marker[which(gyn_nuc_proto_marker$ATG %in% unique_gene_ids_proto),]$gene

gyn.combined_add = AddModuleScore(object = gyn.combined,features = list(unique_gene_ids_proto_gyn),name = 'stress')
DefaultAssay(gyn.combined_add)="RNA"
FeaturePlot(gyn.combined_add,features = 'stress1',split.by ="sample_type",label=F, min.cutoff = -0.1, max.cutoff = 0.5)& scale_colour_gradientn(colours = custom_magma,limits = c(-0.1,0.5))& theme(legend.position = "right")
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da", "#decbe4", "#de77ae",
             "#dfc27d", "#fec44f","#ffed6f", "#ffed6f",
             "#fbb4ae","#bdbdbd")
Idents(gyn.combined_add) <- "sample_type"
VlnPlot(gyn.combined_add,features = 'stress1',pt.size = 0,cols=c("#737373","#addd8e"))
y_max <- max(FetchData(gyn.combined_add, vars = 'stress1',), na.rm = TRUE)
library(ggpubr)

VlnPlot(
  gyn.combined_add,
  features = 'stress1',
  group.by = "sample_type",
  cols = c("#737373", "#addd8e"),
  pt.size = 0,
  y.max = y_max + 0.1 # Add buffer to avoid cutting off p-values
) + 
  stat_compare_means(comparisons = list(c("nuc", "proto")), label = "p.signif")
#gyn.combined_add_nuc=subset(gyn.combined_add,idents = "nuc")
#gyn.combined_add_proto=subset(gyn.combined_add,idents = "proto")
Idents(gyn.combined_add_nuc) <- "new_cell_type"
VlnPlot(gyn.combined_add_nuc,features = 'stress1',pt.size = 0,cols=colour_bar)
Idents(gyn.combined_add_proto) <- "new_cell_type"
VlnPlot(gyn.combined_add_proto,features = 'stress1',pt.size = 0,cols=colour_bar)
library("Nebulosa")
plot_density(gyn.combined_add, "stress1")



filtered_rows_nuc <- gyn_GO_for_plot[gyn_GO_for_plot$Description %in% Descripton_to_keep_nuc, "geneID"]
gene_ids_nuc <- unlist(strsplit(filtered_rows_nuc, "/"))
unique_gene_ids_nuc <- unique(gene_ids_nuc)
unique_gene_ids_nuc_gyn=gyn_nuc_proto_marker[which(gyn_nuc_proto_marker$ATG %in% unique_gene_ids_nuc),]$gene
gyn.combined_add = AddModuleScore(object = gyn.combined_add,features = list(unique_gene_ids_nuc_gyn),name = 'nuc')
DefaultAssay(gyn.combined_add)="RNA"
FeaturePlot(gyn.combined_add,features = 'nuc1',split.by ="sample_type",label=T)& scale_colour_gradientn(colours = custom_magma)& theme(legend.position = "right")
FeaturePlot(gyn.combined_add,features = 'nuc1',split.by ="sample_type",label=F, min.cutoff = 0, max.cutoff = 0.25)&
  scale_colour_gradientn(colours = custom_magma,limits = c(0,0.25))& 
  theme(legend.position = "right")

Idents(gyn.combined_add) <- "new_cell_type"
y_max <- max(FetchData(gyn.combined_add, vars = 'nuc1',), na.rm = TRUE)
VlnPlot(
  gyn.combined_add,
  features = 'nuc1',
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
Idents(gyn.combined_add) <- "sample_type"
gyn.combined_add_nuc=subset(gyn.combined_add,idents = "nuc")
gyn.combined_add_proto=subset(gyn.combined_add,idents = "proto")
Idents(gyn.combined_add_nuc) <- "new_cell_type"
a=VlnPlot(gyn.combined_add_nuc,features = 'nuc1',pt.size = 0,cols=colour_bar)
Idents(gyn.combined_add_proto) <- "new_cell_type"
b=VlnPlot(gyn.combined_add_proto,features = 'nuc1',pt.size = 0,cols=colour_bar)
