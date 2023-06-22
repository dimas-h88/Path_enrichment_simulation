create_res <- function(res, name, ks){
  name <- data.table(pathway = res$pathway,
                                  pval = res$pval,
                                  padj = res$padj,
                                  seed = i,
                                  N = n,
                                  Method = name, 
                                  k = ks)
  return(name)
}