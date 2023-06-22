do_bi_fora <- function(de_table, gene_sets=geneSets){
  nSignGenesFraction <- sum(de_table$padj < 0.05, na.rm = TRUE) / nrow(de_table)
  
  if (nSignGenesFraction > 0.25) {
    genesUp <- rownames(de_table[(de_table$padj) < 0.05 & de_table$stat >= 0, ])
    genesDown <- rownames(de_table[de_table$padj < 0.05 & de_table$stat < 0, ])
    
  } else if (nSignGenesFraction < 0.01) {
    genesUp <- rownames(de_table[de_table$pvalue < 0.05 & de_table$stat >= 0, ])
    genesDown <- rownames(de_table[de_table$pvalue < 0.05 & de_table$stat < 0, ])
  } else {
    genesUp <- rownames(de_table[de_table$padj < 0.05 & de_table$stat >= 0, ])
    genesDown <- rownames(de_table[de_table$padj < 0.05 & de_table$stat < 0, ])
  }
  
  
  
  foraUpRes <- fora(gene_sets, (genesUp), de_table@rownames)
  
  foraDownRes <- fora(gene_sets, (genesDown), de_table@rownames)
  
  # Extract the pathway column as a vector
  #pathwayVec <- as.vector(foraDownRes$pathway)
  
  foraAggrRes <- data.table(pathway=foraUpRes$pathway, pval1=foraUpRes$pval)
  foraAggrRes[, pval2 := foraDownRes[match(foraAggrRes$pathway, foraDownRes$pathway), "pval"]]
  
  aggrPvals <- sapply(seq_len(nrow(foraAggrRes)), function(i){
    return(aggregation::fisher(c(foraAggrRes[i, pval1], foraAggrRes[i, pval2])))
  })
  foraAggrRes[, aggrPval := aggrPvals]
  return(foraAggrRes[, list(pathway, pval=aggrPval)])
}

