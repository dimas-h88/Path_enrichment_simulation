camera_voom_m0 <- function(samp_exprs, samp_meta, dds){
  
  design <- data.frame(samp_meta$group <- rep(c("group1", "group2"), each = n/2))
  
  rownames(design) <- colnames (samp_exprs)
  colnames(design) <- 'Group'
  design_matrix <- model.matrix(~Group, data=design)[,-1]
  design_matrix <- as.matrix(design_matrix)
  
  y <- voom(counts(dds), design = design_matrix, plot = F)
  cameraResults <- camera(y, index = pathways, design = design_matrix, inter.gene.cor=0.01, allow.neg.cor= TRUE)
  
  result_table_camera <- data.table(pathway = rownames(cameraResults),
                                    pval = cameraResults$PValue,
                                    padj = p.adjust(cameraResults$PValue, method = 'BH'),
                                    seed = i,
                                    N = n,
                                    Method = 'CAMERA')
  return(result_table_camera) 
}