fgsea_Multi <- function(res, pathways) {
  rwm <- res$rn
  log_f <- res$log2FoldChange
  names(log_f) <- rwm
  log_f <- log_f[!is.na(log_f)]
  log_f <- log_f[order(log_f)]
  fgseaRes <- fgseaMultilevel(pathways = pathways, 
                              stats    = log_f,
                              minSize = 15,
                              maxSize = 500
  )
  result_table_fgsea <- data.table(pathway = fgseaRes$pathway,
                                   pval = fgseaRes$pval,
                                   padj = fgseaRes$padj,
                                   seed = i,
                                   N = n,
                                   Method = 'FGSEA'
  )
  return (result_table_fgsea)
}