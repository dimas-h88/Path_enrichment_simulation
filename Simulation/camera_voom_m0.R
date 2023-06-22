camera_voom_m0 <- function(counts1, design_matrix, pathways){
  
 
  cameraResults <- camera(counts1, index = pathways, design = design_matrix, inter.gene.cor=0.01, allow.neg.cor= TRUE)
  
  result_table_camera <- data.table(pathway = rownames(cameraResults),
                                    pval = cameraResults$PValue,
                                    padj = p.adjust(cameraResults$PValue, method = 'BH'),
                                    seed = i,
                                    N = n,
                                    Method = 'CAMERA')
  return(result_table_camera) 
}