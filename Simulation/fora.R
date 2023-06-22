fora_res <- function(de_table, gene_sets, ...) {
  
  nSignGenesFraction <- sum(p.adjust(de_table$pval, method = "BH") < 0.05, na.rm = TRUE) / nrow(de_table)
  
  if (nSignGenesFraction > 0.25) {
    genesUp <- de_table[de_table$pval < 0.05 & de_table$stat >= 0, rn]
    genesDown <- de_table[de_table$pval < 0.05 & de_table$stat < 0, rn]
  } else if (nSignGenesFraction < 0.01) {
    genesUp <- de_table[de_table$pval < 0.05 & de_table$stat >= 0, rn]
    genesDown <- de_table[de_table$pval < 0.05 & de_table$stat < 0, rn]
  } else {
    genesUp <- de_table[de_table$pval < 0.05 & de_table$stat >= 0, rn]
    genesDown <- de_table[de_table$pval < 0.05 & de_table$stat < 0, rn]
  }
  
  foraUpRes <- fora(gene_sets, (genesUp), de_table$rn)
  foraDownRes <- fora(gene_sets, genesDown, de_table$rn)
  
  foraAggrRes <- data.table(pathway = foraUpRes$pathway, pval1 = foraUpRes$pval)
  foraAggrRes[, pval2 := foraDownRes[match(foraAggrRes$pathway, foraDownRes$pathway), pval]]
  
  aggrPvals <- sapply(seq_len(nrow(foraAggrRes)), function(i) {
    return(aggregation::fisher(c(foraAggrRes[i, pval1], foraAggrRes[i, pval2])))
  })
  
  foraAggrRes[, aggrPval := aggrPvals]
  return(foraAggrRes[, list(pathway, pval = aggrPval)])
}


do_bi_fora200 <- function(de_table, gene_sets=geneSets, ...){
  genesUp <- de_table[head(order(stat, decreasing = TRUE), 200), rn]
  genesDown <- de_table[tail(order(stat, decreasing = TRUE), 200), rn]
  
  foraUpRes <- fora(gene_sets, genesUp, de_table$rn, ...)
  foraDownRes <- fora(gene_sets, genesDown, de_table$rn, ...)
  
  foraAggrRes <- data.table(pathway=foraUpRes$pathway, pval1=foraUpRes$pval)
  foraAggrRes[, pval2 := foraDownRes[match(foraAggrRes$pathway, foraDownRes$pathway), pval]]
  
  aggrPvals <- sapply(seq_len(nrow(foraAggrRes)), function(i){
    return(aggregation::fisher(c(foraAggrRes[i, pval1], foraAggrRes[i, pval2])))
  })
  
  foraAggrRes[, aggrPval := aggrPvals]
  return(foraAggrRes[, list(pathway, pval=aggrPval)])
  
  
  do_bi_fora <- function(de_table, gene_sets=geneSets, ...){
    nSignGenesFraction <- sum(p.adjust(de_table$pvalue, method = "BH") < 0.05, na.rm = TRUE) / nrow(de_table)
    
    if (nSignGenesFraction > 0.25) {
      genesUp <- de_table[p.adjust(pvalue) < 0.05 & stat >= 0, rn]
      genesDown <- de_table[p.adjust(pvalue) < 0.05 & stat < 0, rn]
    } else if (nSignGenesFraction < 0.01) {
      genesUp <- de_table[pvalue < 0.05 & stat >= 0, rn]
      genesDown <- de_table[pvalue < 0.05 & stat < 0, rn]
    } else {
      genesUp <- de_table[p.adjust(pvalue, method = "BH") < 0.05 & stat >= 0, rn]
      genesDown <- de_table[p.adjust(pvalue, method = "BH") < 0.05 & stat < 0, rn]
    }
    
    foraUpRes <- fora(gene_sets, genesUp, de_table$rn, ...)
    foraDownRes <- fora(gene_sets, genesDown, de_table$rn, ...)
    
    foraAggrRes <- data.table(pathway=foraUpRes$pathway, pval1=foraUpRes$pval)
    foraAggrRes[, pval2 := foraDownRes[match(foraAggrRes$pathway, foraDownRes$pathway), pval]]
    
    aggrPvals <- sapply(seq_len(nrow(foraAggrRes)), function(i){
      return(aggregation::fisher(c(foraAggrRes[i, pval1], foraAggrRes[i, pval2])))
    })
    
    foraAggrRes[, aggrPval := aggrPvals]
    return(foraAggrRes[, list(pathway, pval=aggrPval)])
  }
}
