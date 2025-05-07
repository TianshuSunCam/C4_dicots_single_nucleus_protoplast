library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

Idents(bid.combined) <- "sample_type"
DefaultAssay(bid.combined)="RNA"
options(future.globals.maxSize = 1000 * 1024^2)
all.varied.re.markers <- FindMarkers(bid.combined, ident.1 = "proto", ident.2 = "nuc", verbose = T)
write.table(all.varied.re.markers, 'bidentis_integrated_nuc_vs_proto_marker.tab',col.names = T,row.names = T,quote = F)

###bid
bid_nuc_proto_marker=read.table("./bidentis_integrated_nuc_vs_proto_marker.tab")
###bid_nuc_proto_marker=read.table("./all_markers_re_vaired.tab")
bid_nuc_proto_marker$gene=rownames(bid_nuc_proto_marker)
ara_bid=read.table("./new_names_besthit_bidentis_arabidopsis.out",header = F)
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
write.csv(bid_test_enrich@compareClusterResult,"bid_GO_compareClusterResult_proto_vs_nuc.csv")

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
