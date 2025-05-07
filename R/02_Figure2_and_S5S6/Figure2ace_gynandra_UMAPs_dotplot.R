library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
###gyn
gyn.combined=readRDS("yourpath/gyn.combined_sub.rds")
###cluster annotation and dimplot
Idents(gyn.combined)="new_cell_type"
colour_bar=c("#238443", "#1d91c0",
             "#9970ab","#f1b6da", "#decbe4", "#de77ae",
             "#dfc27d", "#fec44f","#ffed6f", "#ffed6f",
             "#fbb4ae","#bdbdbd")
DimPlot(gyn.combined, reduction = "umap",label=F,group.by = "new_cell_type",cols=colour_bar)
DefaultAssay(gyn.combined)='RNA'
###
#######cluster annotation
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
gyn.combined.core=subset(x = gyn.combined, idents = c("mesophyll","BS","epidermis","guard","phloem","xylem","proliferating"))
saveRDS(gyn.combined, file = "./gyn_combined_core.rds")
gyn.combined=readRDS("./gyn_combined_core.rds")

global_range <- c(-0.5,2.5)

gyn_Dot=DotPlot(object =  gyn.combined.core,
                features = rev(query_gene),scale.by = 'size',group.by = 'new_cell_type',  col.min = -0.5,
                col.max = 2.5,,scale.min = 0, scale.max = 100)+ 
  coord_flip() +theme(axis.text.x = element_text(size = 14,angle = 30,vjust=1,hjust = 1))& scale_colour_gradientn(colours = custom_magma, limits = global_range)
gyn_Dot
#####
Idents(gyn.combined)='new_cell_type'
all.markers <- FindAllMarkers(object = gyn.combined,only.pos = TRUE)
write.csv(all.markers, 'Gg_all_markers_cell_type.csv',col.names = T,row.names = T,quote = F)
Gyn_all_markers=read.csv("./Gg_all_markers_cell_type.csv")
ara_gyn=read.table("./besthit_gyn_to_arabidopsis.out",header = F)
ara_gyn$V1=gsub("\\.mRNA[0-9]+","",ara_gyn$V1)
ara_gyn$V2=gsub("\\.[0-9]+","",ara_gyn$V2)
Gyn_all_markers$ATG=ara_gyn$V2[match(Gyn_all_markers$gene, ara_gyn$V1)]
###library(DOSE)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(tidyverse)
library("org.At.tair.db")
meso_genes=Gyn_all_markers[which(Gyn_all_markers$cluster=='mesophyll'),]$ATG[1:200]
meso_GO = enrichGO(meso_genes, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,readable = TRUE)
head(meso_GO)
dotplot(meso_GO, showCategory=30) + ggtitle("dotplot for mesophyll")

####
BS_genes=Gyn_all_markers[which(Gyn_all_markers$cluster=='BS'),]$ATG[1:200]
BS_GO = enrichGO(BS_genes, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,readable = TRUE)
head(BS_GO)
dotplot(BS_GO, showCategory=30) + ggtitle("dotplot for BS")
###
xylem_genes=Gyn_all_markers[which(Gyn_all_markers$cluster=='xylem'),]$ATG[1:200]
xylem_GO = enrichGO(xylem_genes, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,readable = TRUE)
head(xylem_GO)
dotplot(xylem_GO, showCategory=30) + ggtitle("dotplot for xylem")###


meso_genes=Gyn_all_markers[which(Gyn_all_markers$cluster=='mesophyll'),]$ATG[1:200]
BS_genes=Gyn_all_markers[which(Gyn_all_markers$cluster=='BS'),]$ATG[1:200]
xylem_genes=Gyn_all_markers[which(Gyn_all_markers$cluster=='xylem'),]$ATG[1:200]
phloem_genes=Gyn_all_markers[which(Gyn_all_markers$cluster=='phloem'),]$ATG[1:200]
epidermis_genes=Gyn_all_markers[which(Gyn_all_markers$cluster=='epidermis'),]$ATG[1:200]
guard_genes=Gyn_all_markers[which(Gyn_all_markers$cluster=='guard'),]$ATG[1:200]
proliferating_genes=Gyn_all_markers[which(Gyn_all_markers$cluster=='proliferating'),]$ATG[1:200]
###
query=Gyn_all_markers[which(Gyn_all_markers$cluster=='proliferating'),]$ATG[1:200]
query_GO = enrichGO(query, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP",pAdjustMethod = "BH", pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,readable = TRUE)
head(query_GO)
dotplot(query_GO, showCategory=30) + ggtitle("dotplot for query")
###
enrich_gene_list_for_test <- list('mesophyll'=meso_genes,'BS'=BS_genes,
                                  'xylem'=xylem_genes,'phloem'=phloem_genes,
                                  'epidermis'=epidermis_genes,'guard'=guard_genes,
                                  'proliferating'=proliferating_genes)
gyn_test_enrich <- compareCluster(geneCluster = enrich_gene_list_for_test,ont = "BP", 
                              OrgDb = org.At.tair.db, pAdjustMethod = "BH",keyType="TAIR",
                              pvalueCutoff  = 0.01, qvalueCutoff  = 0.05,
                              fun =  enrichGO)
dotplot(gyn_test_enrich, showCategory=30) + ggtitle("dotplot")+coord_flip()
gyn_GO_for_plot=gyn_test_enrich@compareClusterResult
Descripton_to_keep=c("photosynthesis","photosynthesis, light reaction",
                     "electron transport chain","chlorophyll metabolic process",
                     "NAD(P)H dehydrogenase complex assembly","nonphotochemical quenching",
                     "chloroplast organization","photosynthesis, dark reaction",
                     "glycine metabolic process","serine family amino acid metabolic process",
                     "cellular polysaccharide metabolic process",
                     "cell wall biogenesis","xylem development","response to auxin",
                     "regionalization","pattern specification process","homeostatic process",
                     "monocarboxylic acid biosynthetic process","fatty acid biosynthetic process","cuticle development","wax biosynthetic process",
                     "response to carbon dioxide","stomatal complex development",
                     "cell cycle process","microtubule-based process")

head(gyn_GO_for_plot)
gyn_GO_for_plot=gyn_test_enrich@compareClusterResult
gyn_GO_for_plot=gyn_GO_for_plot[which(gyn_GO_for_plot$Description %in% Descripton_to_keep),]
gyn_GO_for_plot$Description=factor(gyn_GO_for_plot$Description,levels = rev(Descripton_to_keep))

gyn_GO_for_plot$GeneRatio_numeric <- sapply(gyn_GO_for_plot$GeneRatio, function(x) {
  ratio <- as.numeric(unlist(strsplit(x, "/")))  # Split and convert to numeric
  ratio[1] / ratio[2]  # Divide the two numbers
})

gyn_GO=ggplot(gyn_GO_for_plot, aes(x=Cluster, y=Description,color= p.adjust)) +
  geom_point(aes(size=GeneRatio_numeric))+scale_size(limits=c(0,0.5),breaks=c(0.1,0.2,0.3,0.4,0.5))+
  #geom_segment(aes(x = 0, xend = Cluster, y =Description, yend = Description),size = 1,color='grey')+
  scale_colour_gradientn(colours = c( "#440154FF" ,"#414487FF", "#2A788EFF" ,"#22A884FF" ,"#7AD151FF", "#FDE725FF"),limits = c(0,0.0056))+
  theme(text = element_text(size=16)) +
  theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1, vjust =1, face = "plain"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

write.csv(gyn_test_enrich@compareClusterResult,"gyn_GO_compareClusterResult.csv")
