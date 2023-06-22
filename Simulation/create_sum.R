create_sum <- function(name){
  name <- data.table(pathway = character(),
                                            pval = numeric(),
                                            padj = numeric(),
                                            seed = numeric(),
                                            N = numeric(),
                                            Method = character(),
                                            k = numeric())
  
  
  return(name)
}