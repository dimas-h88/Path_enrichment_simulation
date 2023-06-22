diffexpressed_pathways <- function(res, col){
  
  pathways <- res[res[[col]] < 0.05, 'pathway']
  pathways <- (pathways$pathway)
  return(pathways)  
}


report <- function(method, dif_path, dif_path_gained, n, i, pathways, k){
  
  false_positives <- dif_path_gained[!(dif_path_gained %in% dif_path)]
  true_positives <- dif_path_gained[(dif_path_gained %in% dif_path)]
  
  negatives <- pathways[!(names(pathways) %in% dif_path)]
  negatives_gained <- pathways[!(names(pathways) %in% dif_path_gained)]
  
  false_negatives <- negatives_gained[!(negatives_gained %in% negatives)]
  true_negatives <- negatives_gained[(negatives_gained %in% negatives)]
  
  
  FP <- length(false_positives)
  TP <- length(true_positives)
  FN <- length(false_negatives)
  TN <- length(true_negatives)
  
  # if value is Nan, set it to 0
  true_positive_rate <- ifelse(is.nan(TP/(TP+FN)), 0, (TP/(TP+FN)))
  true_negative_rate <- ifelse(is.nan(TN/(FP+TN)), 0, (TN/(FP+TN)))
  false_positive_rate <- ifelse(is.nan(FP/(FP+TN)), 0, (FP/(FP+TN)))
  false_negative_rate <- ifelse(is.nan(FN/(FN+TP)), 0, (FN/(FN+TP)))
  
  false_discovery_rate <- ifelse(is.nan(FP/(FP+TP)), 0, (FP/(FP+TP)))
  
  accuracy <- (FP)/(length(dif_path))
  
  
  dif_path_frac <- length(dif_path_gained)/length(pathways)
  
  report <- data.table(method = method,
                       true_positive_rate = true_positive_rate,
                       true_negative_rate = true_negative_rate,
                       false_positive_rate = false_positive_rate,
                       false_negative_rate = false_negative_rate,
                       false_discovery_rate = false_discovery_rate,
                       accuracy = accuracy,
                       N = n,
                       seed = i,
                       k = k,
                       Fraction_of_differentially_expressed_pathways = dif_path_frac)
  
  
  return(report)
  
}