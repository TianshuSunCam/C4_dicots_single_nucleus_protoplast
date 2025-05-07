library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
gyn.combined=readRDS("./gyn.combined_sub.rds")

Idents(gyn.combined) <- "sample_type"
DefaultAssay(gyn.combined)="RNA"
options(future.globals.maxSize = 1000 * 1024^2)
all.varied.re.markers <- FindMarkers(gyn.combined, ident.1 = "proto", ident.2 = "nuc", verbose = T)
write.table(all.varied.re.markers, 'gynentis_integrated_nuc_vs_proto_marker.tab',col.names = T,row.names = T,quote = F)

###gyn
gyn_nuc_proto_marker=read.table("./gynandra_integrated_nuc_vs_proto_marker.tab")
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
proto_genes_GO = enrichGO(proto_genes, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,readable = TRUE)
head(proto_genes_GO)
dotplot(proto_genes_GO, showCategory=30) + ggtitle("dotplot for protoplast")

nuc_genes=gyn_nuc_proto_marker[which(gyn_nuc_proto_marker$avg_log2FC<0),]$ATG
nuc_genes_GO = enrichGO(nuc_genes, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,readable = TRUE)
head(nuc_genes_GO)
dotplot(nuc_genes_GO, showCategory=30) + ggtitle("dotplot for nuclei")
###
enrich_gene_list_for_test <- list('protoplast'=proto_genes,'nuclei'=nuc_genes)
gyn_test_enrich <- compareCluster(geneCluster = enrich_gene_list_for_test,ont = "BP", 
                                  OrgDb = org.At.tair.db, pAdjustMethod = "BH",keyType="TAIR",
                                  pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                                  fun =  enrichGO)
dotplot(gyn_test_enrich, showCategory=30) + ggtitle("dotplot")+coord_flip()
gyn_GO_for_plot=gyn_test_enrich@compareClusterResult
write.csv(gyn_test_enrich@compareClusterResult,"gyn_GO_compareClusterResult_proto_vs_nuc.csv")

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



Descripton_to_keep=c(Descripton_to_keep_nuc,rev(Descripton_to_keep_proto))

gyn_GO_for_plot=gyn_test_enrich@compareClusterResult
gyn_GO_for_plot=gyn_GO_for_plot[which(gyn_GO_for_plot$Description %in% Descripton_to_keep),]
gyn_GO_for_plot$Description=factor(gyn_GO_for_plot$Description,levels = c(Descripton_to_keep))

gyn_GO_for_plot$GeneRatio_numeric <- sapply(gyn_GO_for_plot$GeneRatio, function(x) {
  ratio <- as.numeric(unlist(strsplit(x, "/")))  # Split and convert to numeric
  ratio[1] / ratio[2]  # Divide the two numbers
})
ggplot(gyn_GO_for_plot, aes(x=Description, y=Cluster,color= p.adjust)) +
  geom_point(aes(size=GeneRatio_numeric))+
#scale_size(limits=c(0,0.06),breaks=c(0.01,0.02,0.03,0.04,0.05,0.06))+
  #geom_segment(aes(x = Description, xend = Description, y = 0, yend = Cluster),size = 1,color='#d9d9d9')+
  scale_colour_gradientn(colours = c( "#440154FF" ,"#414487FF", "#2A788EFF" ,"#22A884FF" ,"#7AD151FF", "#FDE725FF"))+
  theme(text = element_text(size=16)) +theme(legend.position="top")+
  theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1, vjust =1, face = "plain"),
        panel.grid.major = element_line(color = "#969696", size = 0.5,linetype = 2), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

