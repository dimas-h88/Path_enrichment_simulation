fgsea_voom_m0 <- function(samp_exprs, samp_meta, dds){
  design <- data.frame(samp_meta$Group <- rep(c("group1", "group2"), each = n/2))
  
  rownames(design) <- colnames (exprs)
  colnames(design) <- 'Group'
  design_matrix <- model.matrix(~Group, data=design)[,-1]
  design_matrix <- as.matrix(design_matrix)
  
  y <- voom(counts(dds), design = design_matrix, plot = F)
  fit <- lmFit(y, design_matrix)
  fit <- eBayes(fit)
  top.table <- topTable(fit, sort.by = "P", n = Inf)
  t <- top.table$t
  names(t) <- rownames(top.table)
  t <- t[order(t)]
  fgsea_voom <- fgseaMultilevel(pathways = pathways, 
                                stats    = t,
                                minSize = 15,
                                maxSize = 500)
  
  result_table_f_voom <- data.table(pathway = fgsea_voom$pathway,
                                    pval = fgsea_voom$pval,
                                    padj = fgsea_voom$padj,
                                    seed = i,
                                    N = n,
                                    Method = 'FGSEA_voom')
return(result_table_f_voom)
}