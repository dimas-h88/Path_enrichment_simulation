fgsea_Multi <- function(res, pathways) {
  rwm <- res@rownames
  stats <- res$log2FoldChange
  names(stats) <- rwm
 # stats <- stats[!is.na(stats)]
#  stats <- stats[order(stats)]
  fgseaRes <- fgseaMultilevel(pathways = pathways, 
                              stats    = stats
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