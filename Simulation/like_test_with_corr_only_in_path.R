library(limma)
library(GEOquery)
library(edgeR)
library(phantasus)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(data.table)
library(reactome.db)
library(fgsea)
library(aggregation)
library(data.table)
library(stats)
library(msigdbr)
library(collapse)
library(insight)
library(data.table)
library(Matrix)
library(mvnfast)
source('do_bi_fora.R')
source('fgsea_Multi.R')
source('camera_voom_m0.R')
source('fgsea_voom_m0.R')
source('fgseaResample.R')
source('create_sum.R')
source('creare_res.R')
source('report.R')
source('fora.R')

#BiocParallel::register(BiocParallel::SerialParam())

summary_report <- data.table(method = character(),
                             true_positive_rate = numeric(),
                             true_negative_rate = numeric(),
                             false_positive_rate = numeric(),
                             false_negative_rate = numeric(),
                             false_discovery_rate = numeric(),
                             accuracy = numeric(),
                             N = numeric(),
                             seed = numeric(),
                             k = numeric(),
                             Fraction_of_differentially_expressed_pathways = numeric())


#Function which creates summary tables (stored at 'create_sum.R')
summary_table_SE <- create_sum('summary_table_SE')
summary_table_camera <- create_sum('summary_table_camera')
summary_table_fora <- create_sum('summary_table_fora')
summary_table_fgsea <- create_sum('summary_table_fgsea')

# Importing raw data to take real names
raw_data <- read.gct("New_no_sym.gct")
# Taking summary count matrix
counts <-raw_data@assayData$exprs
entrez_ids <- select(org.Hs.eg.db, keys=rownames(counts), columns="ENTREZID", keytype="ENSEMBL",  uniqueRows=TRUE)
entrez_ids <- entrez_ids[!duplicated(entrez_ids$ENSEMBL),] #deleting multiple conversion
na_rows <- is.na(entrez_ids$ENTREZID) # deleting NA names
counts <- counts[!na_rows,]# deleting rows with NA ENTREZIDs
entrez_ids <- entrez_ids[!na_rows,]
rownames(counts) <- entrez_ids$ENTREZID# Changing names for further analysis
#Filtering low expressed genes (not nessesarry in simulation case)
counts <- counts[complete.cases(rownames(counts)),]
counts <- counts[order(-rowSums(counts)),]
counts <- head(counts, 1000)
gene_names <- as.character(rownames(counts))
n_genes <- length(gene_names)
gene_names <- as.character(gene_names[!duplicated(gene_names)])
#Importing pathways, containing genes from our real world matrix
msigdbr_df <- msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'GO:BP')
msigdbr_list = split(x = msigdbr_df$entrez_gene, f = msigdbr_df$gs_name)

# Filter the list to keep only the elements that have intersections with entrez_ids,
## taking only unique values(nessessarry for GO:pathways, containing duplicates in 1 pathway)
pathways <- lapply(msigdbr_list, intersect, y = gene_names)
pathways <- pathways[!is.na(pathways)]
pathways <- lapply(pathways, unique)
pathways <- lapply(pathways, as.character)

#Filtering pathways by length
pathways <- pathways[lapply(pathways, length) >= 15]
pathways <- pathways[lapply(pathways, length) < 500]



pathway_len <- as.vector(lapply(pathways, length))
path_len <- vector()
for (v in 1:length(pathways)) {
  path_len[v] <- pathway_len[[v]]
}
generated_pathways <- list()
for (v in 1:length(pathways)){
  set.seed(v)
  genes_in_pathway <- as.character(sample(gene_names, path_len[v]))
  generated_pathways[[v]] <- (genes_in_pathway)
}

#char_vector <- as.character(seq(from = 1, to = length(generated_pathways)))
#char_vector <- paste("P", char_vector, sep = "")
generated_pathways <- setNames(generated_pathways, names(pathways))
pathways <- generated_pathways




simulate_vector <- function(n_genes, gene_names, 
                            noize_sd, seed, log2FC_sc = NULL, 
                            diffexp_genes = NULL) 
{
  expression_vector <- rep(10, n_genes)
  
  noize <- rnorm(n_genes, 0, noize_sd)
  expression_vector <- expression_vector + noize
  
  if (!is.null(diffexp_genes) & !is.null(log2FC_sc)) {
    log2FC <- rep(0, n_genes)
    names(log2FC) <- gene_names
    
    # Taking genes from exact pathways and changing log2FC of them from 0 to 1
    log2FC[diffexp_genes] <- log2FC_sc
    expression_vector <- as.numeric(expression_vector + log2FC)
  }
  
  return(expression_vector)
}

#plotEnrichment(pathways[["GOBP_ACTIN_FILAMENT_ORGANIZATION"]], stats)

#sigma <- matrix(rnorm(n_genes^2, 0, 0.01), nrow = n_genes)

for (k in c(0.02, 0.05, 0.1,  0.5, 1, 2, 10)){
  for (n in c(10)) {
    for (i in 1:10) {
      
      
      
      set.seed(i)
      diffexp_genes <- as.character(pathways[["GOBP_ACTIN_CYTOSKELETON_REORGANIZATION"]])
      samp_path <- sample(pathways, 10)
      samp_genes <- unique(unlist(samp_path))
      samp_genes <- unique(c(samp_genes, diffexp_genes))
      #correlated_pathways_genes <- (unique(c(pathways[["GOBP_ACTIN_FILAMENT_ORGANIZATION"]], pathways[["GOBP_CELLULAR_RESPONSE_TO_STEROID_HORMONE_STIMULUS"]], pathways[["GOBP_CLATHRIN_DEPENDENT_ENDOCYTOSIS"]] )))
    
      
      
      expected_expression <- rep(10, n_genes)
      names(expected_expression) <- gene_names
      expected_expression_2 <- expected_expression
      expected_expression_2[diffexp_genes] <- 11
      
      sigma <- matrix(runif(n_genes^2, 0, 0.01), nrow = n_genes)
      colnames(sigma) <- as.character(gene_names)
      rownames(sigma) <- as.character(gene_names)
      
      for(m in lapply(samp_path, unlist)){
        for (z in m)
        {
          for (j in m)
          {
            sigma[as.character(z), as.character(j)] <- k
          }
        }
        
        for (z in m)
        {
          sigma[z, z] <- k/2
        }
      }
      
      sigma <- Matrix::nearPD(sigma)$mat
      
      set.seed(i)
      
      #sigma <- (sigma + t(sigma)) / 2
      #sigma <- nearPD(sigma)$mat
      
      
      
      
      simulate_corr_vector <- function(n_genes, 
                                       
                                       sigma,
                                       expression_mean_vector,
                                       log2FC_n = NULL, 
                                       diffexp_genes = NULL)
      {
        # Generate random expression vector
        expression_vector <- mvnfast::rmvn(1, expression_mean_vector, as.matrix(sigma))
        
        if (!is.null(diffexp_genes) & !is.null(log2FC_n)) {
          log2FC <- rep(0, n_genes)
          names(log2FC) <- gene_names
          
          # Taking genes from exact pathways and changing log2FC of them from 0 to 1
          log2FC[diffexp_genes] <- log2FC_n
          expression_vector <- expression_vector + log2FC
        }  
        return(expression_vector)
      }
      
      
      
      
      
      
      pathways <- generated_pathways
      pathways <- lapply(pathways, as.character)
      
      print(i)
      print(n)
      
      
      
      # Simulating clear vector of expressions 
      matrix_expression_state_1 <- matrix(nrow = n_genes, ncol = n)
      rownames(matrix_expression_state_1) <- gene_names
      
      # simulate samples and append to the matrix
      for (j in 1:n) {
        set.seed(j)
        expression_state_1 <- simulate_corr_vector(n_genes,
                                                   sigma = sigma,
                                                   expected_expression
        )
        matrix_expression_state_1[, j] <- expression_state_1
      }
      
      
      matrix_expression_state_2 <- matrix(nrow = length(gene_names), ncol = n)
      rownames(matrix_expression_state_2) <- gene_names
      for (j in 1:n) {
        set.seed(j+n)
        expression_state_2 <- simulate_corr_vector(n_genes,
                                                   sigma = sigma,
                                                   expected_expression_2
        )
        matrix_expression_state_2[, j] <- expression_state_2
      }
      
      set.seed(i)
      
      
      colnames(matrix_expression_state_1) <- rep('State_1', n) 
      colnames(matrix_expression_state_2) <- rep('State_2', n) 
      count_matrix <- (cbind(matrix_expression_state_1, matrix_expression_state_2))
      count_matrix <- normalizeBetweenArrays(count_matrix)
      # Extracting pathways, which containing true differentially expressed genes(Which
      # log2FC was changed)
      
      pathways_containing_diffexp_genes <- names(pathways[lengths(lapply(pathways, intersect, diffexp_genes))])
      
      # STARTING Limma
      design <- data.frame(Group = factor(c(rep("group1", n), rep("group2", n))))
      design_matrix <- (model.matrix(~ Group, data = design))
      
      fit <- lmFit(count_matrix, design_matrix)
      fit2 <- eBayes(fit)
      de <- topTable(fit2, adjust.method = "BH", number=Inf, confint = (pnorm(1)-0.5)*2)
      de <- data.table(as.data.frame(de), keep.rownames = TRUE)
      de <- de[, list(rn, 
                      log2FoldChange=logFC, lfcSE=(CI.R-CI.L)/2, 
                      pval=P.Value, stat=t, padj = adj.P.Val)]
      
      # STARTING FGSEA
      
      stats <- de[, setNames(log2FoldChange, rn)]
      
      fgseaSERes <- fgseaMultilevelSE(pathways = pathways, 
                                      stats, se=de$lfcSE)
      
      # Creating result table for this N and seed and adding it to summary table 
      result_table_fgsea_SE <- create_res(fgseaSERes, 'FGSEA_SE',  ks = k)
      summary_table_SE <- rbind(summary_table_SE, result_table_fgsea_SE)
      
      # Taking pathways, with pval < 0.05 after FGSEA_SE (function stored in report.R)
      diffexpressed_pathways_fgsea_SE <- diffexpressed_pathways(result_table_fgsea_SE, 'padj')
      
      # Calculating FP, TN, FN, TP and accuracy and adding it to summary report
      report_fgsea_SE <- report('fgsea_SE', pathways_containing_diffexp_genes, 
                                diffexpressed_pathways_fgsea_SE, n, i, pathways, k)
      summary_report <- rbind(report_fgsea_SE, summary_report)
      
      
      
      # STARTING CAMERA
      
      cameraRES <- camera_voom_m0(count_matrix, design_matrix, pathways)
      
      result_table_camera <- create_res(cameraRES, 'camera', ks = k)
      summary_table_camera <- rbind(summary_table_camera, result_table_camera)
      
      
      diffexpressed_pathways_camera <- diffexpressed_pathways(result_table_camera, 'padj')
      
      report_camera <- report('CAMERA', pathways_containing_diffexp_genes,
                              diffexpressed_pathways_camera, n, i, pathways, k)
      summary_report <- rbind(report_camera, summary_report)
      
      
      # STARTING FORA
      
      
      foraRes <- fora_res(de, pathways)
      
      
      result_table_fora <- data.table(pathway = foraRes$pathway,
                                      pval = foraRes$pval,
                                      padj = p.adjust(foraRes$pval, method = 'BH'),
                                      seed = i,
                                      N = n,
                                      k = k,
                                      Method = 'fora')
      summary_table_fora <- rbind(summary_table_fora, result_table_fora)
      
      
      diffexpressed_pathways_fora <- diffexpressed_pathways(result_table_fora, 'padj')
      
      report_fora <- report('fora', pathways_containing_diffexp_genes, diffexpressed_pathways_fora,
                            n, i, pathways, k)
      summary_report <- rbind(report_fora, summary_report)
      
      # STARTING FGSEA MULTILEVEL  
      
      fstats <- de[, setNames(stat, rn)]
      fgseaRes <- fgseaMultilevel(pathways = pathways, 
                                  fstats)
      
      
      
      # Creating result table http://127.0.0.1:10357/graphics/plot_zoom_png?width=1004&height=684for this N and seed and adding it to summary table 
      result_table_fgsea <- create_res(fgseaRes, 'FGSEA',  ks = k)
      summary_table_fgsea <- rbind(summary_table_fgsea, result_table_fgsea)
      
      # Taking pathways, with pval < 0.05 after FGSEA_SE (function stored in report.R)
      diffexpressed_pathways_fgsea <- diffexpressed_pathways(result_table_fgsea, 'padj')
      
      # Calculating FP, TN, FN, TP and accuracy and adding it to summary report
      report_fgsea <- report('fgsea', pathways_containing_diffexp_genes,
                             diffexpressed_pathways_fgsea, n, i, pathways, k)
      summary_report <- rbind(report_fgsea, summary_report)
      
    }}}
write.csv(summary_table_camera, 'summary_table_camera_padj.csv')
write.csv(summary_table_fgsea, 'summary_table_fgsea_padj.csv')
write.csv(summary_table_fora, 'summary_table_fora_padj.csv')
write.csv(summary_table_SE, 'summary_table_SE_padj.csv')

write.csv(summary_report, 'camera_like_report_padj.csv')


