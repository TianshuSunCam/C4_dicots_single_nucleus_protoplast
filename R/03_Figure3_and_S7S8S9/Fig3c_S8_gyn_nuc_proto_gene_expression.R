library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
###gyn
gyn.combined=readRDS("./gyn.combined_sub.rds")
Idents(gyn.combined)='sample_type'
gyn.combined.nuc=subset(x = gyn.combined, idents = "nuc")
saveRDS(gyn.combined.nuc, file = "./gyn.combined.nuc.rds")
gyn.combined.proto=subset(x = gyn.combined, idents = "proto")
saveRDS(gyn.combined.proto, file = "./gyn.combined.proto.rds")

###read file
gyn.combined.proto=readRDS("./gyn.combined.proto.rds")
gyn.combined.nuc=readRDS("./gyn.combined.nuc.rds")
####
nuc_gyn_counts <- GetAssayData(gyn.combined.nuc,  assay="RNA")   
nuc_gyn_genes.percent.expression <- rowMeans(nuc_gyn_counts>0 )*100   
nuc_gyn_genes.percent.expression
nuc_gyn_genes.filter <- names(nuc_gyn_genes.percent.expression[nuc_gyn_genes.percent.expression>0])
nuc_gyn_genes.filter_1p <- names(nuc_gyn_genes.percent.expression[nuc_gyn_genes.percent.expression>1])  #select genes expressed in at least 1% of cells
nuc_gyn_genes.filter_3p <- names(nuc_gyn_genes.percent.expression[nuc_gyn_genes.percent.expression>3])  #select genes expressed in at least 3% of cells
length(nuc_gyn_genes.filter_1p)
#nuc_gyn_counts.sub <- counts[nuc_gyn_genes.filter,]
#new_seurat_object <- CreateSeuratObject(counts=nuc_gyn_)
nuc_gyn_avg_nCount=mean(gyn.combined.nuc$nCount_RNA)
nuc_gyn_avg_nFeature=mean(gyn.combined.nuc$nFeature_RNA)
######
proto_gyn_counts <- GetAssayData(gyn.combined.proto,  assay="RNA")   
proto_gyn_genes.percent.expression <- rowMeans(proto_gyn_counts>0 )*100   
proto_gyn_genes.percent.expression
proto_gyn_genes.filter <- names(proto_gyn_genes.percent.expression[proto_gyn_genes.percent.expression>0])
proto_gyn_genes.filter_1p <- names(proto_gyn_genes.percent.expression[proto_gyn_genes.percent.expression>1])  #select genes expressed in at least 1% of cells
proto_gyn_genes.filter_3p <- names(proto_gyn_genes.percent.expression[proto_gyn_genes.percent.expression>3])  #select genes expressed in at least 3% of cells
length(proto_gyn_genes.filter)
#proto_gyn_counts.sub <- counts[proto_gyn_genes.filter,]
#new_seurat_object <- CreateSeuratObject(counts=proto_gyn_)
proto_gyn_avg_nCount=mean(gyn.combined.proto$nCount_RNA)
proto_gyn_avg_nFeature=mean(gyn.combined.proto$nFeature_RNA)
######
barplot_gyn=data.frame(type=c(rep("nuclei",2),rep("protoplast",2)),
                       category=c("avg_nFeature","all_detected",
                                  "avg_nFeature","all_detected"),
                       number=c(nuc_gyn_avg_nFeature,
                                length(nuc_gyn_genes.filter),proto_gyn_avg_nFeature,
                                length(proto_gyn_genes.filter)))
barplot_gyn$type=factor(barplot_gyn$type,levels = c("protoplast","nuclei"))
barplot_gyn$category=factor(barplot_gyn$category,levels = c("avg_nFeature","all_detected"))

ggplot(barplot_gyn,aes(x=number,y=type,fill=type))+
  geom_bar(stat='identity',color='gray')+
  geom_text(aes(label=round(number,0)), size=4, position=position_stack(vjust=.5))+
  facet_wrap(category~.,scales = 'free_x', ncol = 1)+
  scale_fill_manual(values=c("#addd8e","#bdbdbd"))+
  #scale_y_continuous(labels = scales::percent_format())+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
#######

#######
#####
###gene_number

gyn_average_expression_proto_unfilter=data.frame(expression=rowMeans(gyn.combined.proto[['RNA']]@data))
gyn_average_expression_proto_unfilter$gene=rownames(gyn_average_expression_proto_unfilter)
gyn_average_expression_proto=data.frame(gyn_average_expression_proto_unfilter[proto_gyn_genes.filter,])
gyn_average_expression_proto$color_nuc <- ifelse(row.names(gyn_average_expression_proto) %in% nuc_gyn_genes.filter, "in", "out")
gyn_average_expression_proto$ranking <- rank(gyn_average_expression_proto$expression)
head(gyn_average_expression_proto)

ggplot(gyn_average_expression_proto, aes(x = ranking, y =expression))+geom_point(aes(color = color_nuc), size =1) + scale_colour_manual(values = c("#bdbdbd","#e31a1c"))+theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank()) 


ggplot(gyn_average_expression_proto, aes(x=ranking, fill=color_nuc)) +
  geom_density(alpha=0.4) +scale_fill_manual(values=c("#addd8e","#737373"))+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

####
mean(gyn_average_expression_proto[which(gyn_average_expression_proto$color_nuc=='out'),1])
mean(gyn_average_expression_proto[which(gyn_average_expression_proto$color_nuc=='in'),1])
gyn_average_expression_proto$color_nuc=factor(gyn_average_expression_proto$color_nuc,levels = c("out","in"))

ggplot(gyn_average_expression_proto, aes(x=color_nuc, y=log2(expression))) +
  geom_violin(aes(fill=color_nuc),alpha=1)+scale_fill_manual(values = c("#bdbdbd","#addd8e"))+
  geom_boxplot(width=0.1,fill="white",alpha=1)+theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
#####
gyn_number_in_proto=gyn_average_expression_proto$expression[ gyn_average_expression_proto$color_nuc== "in"]
gyn_number_out_proto=gyn_average_expression_proto$expression[ gyn_average_expression_proto$color_nuc== "out"]
ttest_result_proto <- t.test(gyn_number_in_proto,gyn_number_out_proto)
ttest_result_proto
###marker gene
Idents(gyn.combined) <- "sample_type"
DefaultAssay(gyn.combined)="RNA"
options(future.globals.maxSize = 1000 * 1024^2)
all.varied.re.markers <- FindMarkers(gyn.combined, ident.1 = "proto", ident.2 = "nuc", verbose = T)
write.csv(all.varied.re.markers, 'gynentis_integrated_nuc_vs_proto_marker.csv',col.names = T,row.names = T,quote = F)