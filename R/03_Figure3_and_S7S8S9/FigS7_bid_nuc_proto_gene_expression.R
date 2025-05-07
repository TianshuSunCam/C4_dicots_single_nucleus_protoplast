library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
###bid
bid.combined=readRDS("./bidentis_combined_integrated_05.rds")
Idents(bid.combined)='sample_type'
bid.combined.nuc=subset(x = bid.combined, idents = "nuc")
bid.combined.proto=subset(x = bid.combined, idents = "proto")
saveRDS(bid.combined.nuc, file = "./bid.combined.nuc.rds")
saveRDS(bid.combined.proto, file = "./bid.combined.proto.rds")

###read file
bid.combined.proto=readRDS("bid.combined.proto.rds")
bid.combined.nuc=readRDS("bid.combined.nuc.rds")
####
nuc_bid_counts <- GetAssayData(bid.combined.nuc,  assay="RNA")   
nuc_bid_genes.percent.expression <- rowMeans(nuc_bid_counts>0 )*100   
nuc_bid_genes.percent.expression
nuc_bid_genes.filter <- names(nuc_bid_genes.percent.expression[nuc_bid_genes.percent.expression>0])
nuc_bid_genes.filter_1p <- names(nuc_bid_genes.percent.expression[nuc_bid_genes.percent.expression>1])  #select genes expressed in at least 1% of cells
nuc_bid_genes.filter_3p <- names(nuc_bid_genes.percent.expression[nuc_bid_genes.percent.expression>3])  #select genes expressed in at least 3% of cells
length(nuc_bid_genes.filter_1p)
#nuc_bid_counts.sub <- counts[nuc_bid_genes.filter,]
#new_seurat_object <- CreateSeuratObject(counts=nuc_bid_)
nuc_bid_avg_nCount=mean(bid.combined.nuc$nCount_RNA)
nuc_bid_avg_nFeature=mean(bid.combined.nuc$nFeature_RNA)
######
proto_bid_counts <- GetAssayData(bid.combined.proto,  assay="RNA")   
proto_bid_genes.percent.expression <- rowMeans(proto_bid_counts>0 )*100   
proto_bid_genes.percent.expression
proto_bid_genes.filter <- names(proto_bid_genes.percent.expression[proto_bid_genes.percent.expression>0])
proto_bid_genes.filter_1p <- names(proto_bid_genes.percent.expression[proto_bid_genes.percent.expression>1])  #select genes expressed in at least 1% of cells
proto_bid_genes.filter_3p <- names(proto_bid_genes.percent.expression[proto_bid_genes.percent.expression>3])  #select genes expressed in at least 3% of cells
length(proto_bid_genes.filter)
#proto_bid_counts.sub <- counts[proto_bid_genes.filter,]
#new_seurat_object <- CreateSeuratObject(counts=proto_bid_)
proto_bid_avg_nCount=mean(bid.combined.proto$nCount_RNA)
proto_bid_avg_nFeature=mean(bid.combined.proto$nFeature_RNA)
######
barplot_bid=data.frame(type=c(rep("nuclei",2),rep("protoplast",2)),
                       category=c("avg_nFeature","all_detected",
                                  "avg_nFeature","all_detected"),
                       number=c(nuc_bid_avg_nFeature,
                                length(nuc_bid_genes.filter),proto_bid_avg_nFeature,
                                length(proto_bid_genes.filter)))
barplot_bid$type=factor(barplot_bid$type,levels = c("protoplast","nuclei"))
barplot_bid$category=factor(barplot_bid$category,levels = c("avg_nFeature","all_detected"))

ggplot(barplot_bid,aes(x=number,y=type,fill=type))+
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

bid_average_expression_proto_unfilter=data.frame(expression=rowMeans(bid.combined.proto[['RNA']]@data))
bid_average_expression_proto_unfilter$gene=rownames(bid_average_expression_proto_unfilter)
bid_average_expression_proto=data.frame(bid_average_expression_proto_unfilter[proto_bid_genes.filter,])
bid_average_expression_proto$color_nuc <- ifelse(row.names(bid_average_expression_proto) %in% nuc_bid_genes.filter, "in", "out")
bid_average_expression_proto$ranking <- rank(bid_average_expression_proto$expression)
head(bid_average_expression_proto)

ggplot(bid_average_expression_proto, aes(x = ranking, y =expression))+geom_point(aes(color = color_nuc), size =1) + scale_colour_manual(values = c("#bdbdbd","#e31a1c"))+theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank()) 


ggplot(bid_average_expression_proto, aes(x=ranking, fill=color_nuc)) +
  geom_density(alpha=0.4) +scale_fill_manual(values=c("#addd8e","#737373"))+
  theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

####
mean(bid_average_expression_proto[which(bid_average_expression_proto$color_nuc=='out'),1])
mean(bid_average_expression_proto[which(bid_average_expression_proto$color_nuc=='in'),1])
bid_average_expression_proto$color_nuc=factor(bid_average_expression_proto$color_nuc,levels = c("out","in"))

ggplot(bid_average_expression_proto, aes(x=color_nuc, y=log2(expression))) +
  geom_violin(aes(fill=color_nuc),alpha=1)+scale_fill_manual(values = c("#bdbdbd","#addd8e"))+
  geom_boxplot(width=0.1,fill="white",alpha=1)+theme_bw() +
  theme(axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


bid_number_in_proto=bid_average_expression_proto$expression[ bid_average_expression_proto$color_nuc== "in"]
bid_number_out_proto=bid_average_expression_proto$expression[ bid_average_expression_proto$color_nuc== "out"]
ttest_result_proto <- t.test(bid_number_in_proto,bid_number_out_proto)
ttest_result_proto
###marker gene
Idents(bid.combined) <- "sample_type"
DefaultAssay(bid.combined)="RNA"
options(future.globals.maxSize = 1000 * 1024^2)
all.varied.re.markers <- FindMarkers(bid.combined, ident.1 = "proto", ident.2 = "nuc", verbose = T)
write.csv(all.varied.re.markers, 'bidentis_integrated_nuc_vs_proto_marker.csv',col.names = T,row.names = T,quote = F)