C1 <- object$score[which(object$cluster == '1')]
C2 <- object$score[which(object$cluster == '2')]
wilcox.test(C1, C2)


vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(gyn.combined_add, features = signature,
            group.by = "sample_type", cols=c("#737373","#addd8e"),pt.size = 0,
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 14, height = 8)
}

gene_sig <- c("stress1")
comparisons <- list(c("nuc", "proto"))
vp_case1(gene_signature = gene_sig, file_name = "./gene_sig_2", test_sign = comparisons)
