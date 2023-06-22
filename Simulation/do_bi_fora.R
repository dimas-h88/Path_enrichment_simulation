do_bi_fora <- function(de_table, gene_sets=geneSets){
  
  de_table <- de
  gene_sets <- pathways
  nSignGenesFraction <- sum(de_table$padj < 0.05, na.rm = TRUE) / nrow(de_table)
  
  if (nSignGenesFraction > 0.25) {
    genesUp <- (de_table[(de_table$padj) < 0.05 & de_table$stat >= 0, ])$rn
    genesDown <- (de_table[de_table$padj < 0.05 & de_table$stat < 0, ])$rn
    print('if')
    
  } else if (nSignGenesFraction < 0.01) {
    genesUp <- (de_table[de_table$pval < 0.05 & de_table$stat >= 0, ])$rn
    genesDown <- (de_table[de_table$pval < 0.05 & de_table$stat < 0, ])$rn
    print('else if')
  } else {
    genesUp <- (de_table[de_table$padj < 0.05 & de_table$stat >= 0, ])$rn
    genesDown <- (de_table[de_table$padj < 0.05 & de_table$stat < 0, ])$rn
    print('else')
  }
  
  
  
  foraUpRes <- fora(gene_sets, (genesUp), de_table$rn)
  
  foraDownRes <- fora(gene_sets, (genesDown), de_table$rn)
  
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


do_bi_fora_1 <- function(de_table, gene_sets=geneSets){
  de_table <- de
  gene_sets <- pathways
  nSignGenesFraction <- sum(de_table$padj < 0.05, na.rm = TRUE) / nrow(de_table)
  
  if (nSignGenesFraction > 0.25) {
    genesUp <- rownames(de_table[(de_table$padj) < 0.05 & de_table$stat >= 0, ])
    genesDown <- rownames(de_table[de_table$padj < 0.05 & de_table$stat < 0, ])
    print('if')
    
  } else if (nSignGenesFraction < 0.01) {
    genesUp <- rownames(de_table[de_table$pval < 0.05 & de_table$stat >= 0, ])
    genesDown <- rownames(de_table[de_table$pval < 0.05 & de_table$stat < 0, ])
    print('else if')
  } else {
    genesUp <- rownames(de_table[de_table$padj < 0.05 & de_table$stat >= 0, ])
    genesDown <- rownames(de_table[de_table$padj < 0.05 & de_table$stat < 0, ])
    print('else')
  }
  
  
  
  foraUpRes <- fora(gene_sets, (genesUp), de_table$rn)
  
  foraDownRes <- fora(gene_sets, (genesDown), de_table$rn)
  
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

