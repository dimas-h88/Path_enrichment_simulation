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
source('do_bi_fora.R')
source('fgsea_Multi.R')
source('camera_voom_m0.R')
source('fgsea_voom_m0.R')
source('fgseaResample.R')
source('create_sum.R')
source('creare_res.R')
source('report.R')
source('fora.R')

summary_report <- data.table(method = character(),
                             true_positive_rate = numeric(),
                             true_negative_rate = numeric(),
                             false_positive_rate = numeric(),
                             false_negative_rate = numeric(),
                             false_discovery_rate = numeric(),
                             accuracy = numeric(),
                             N = numeric(),
                             seed = numeric(),
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
gene_names <- rownames(counts)
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
for (i in 1:length(pathways)) {
  path_len[i] <- pathway_len[[i]]
}
generated_pathways <- list()
for (i in 1:length(pathways)){
  set.seed(i)
  genes_in_pathway <- as.character(sample(gene_names, path_len[i]))
  generated_pathways[[i]] <- (genes_in_pathway)
}

char_vector <- as.character(seq(from = 1, to = length(generated_pathways)))
#char_vector <- paste("P", char_vector, sep = "")
generated_pathways <- setNames(generated_pathways, names(pathways))


simulate_vector <- function(n_genes, gene_names, 
                            noize_sd, seed, log2FC_sc = NULL, 
                            diffexp_genes = NULL) 
{
  expression_vector <- rep(10, n_genes)
  
  set.seed(seed)
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



set.seed(1)
sigma <- matrix(runif(n_genes^2, 0, 0.5), nrow = n_genes)
#sigma <- matrix(rnorm(n_genes^2, mean = 0, sd = 0.25), nrow = n_genes)
#sigma <- as.matrix(sigma)
sigma <- nearPD(sigma)$mat

for (i in 1:n_genes)
  {
   sigma[i, i] <- 0.5
}
#sigma <- nearPD(sigma)$mat

# Every sym matrix is positive
# Add dispersion (diagonal) +- 0.25
# Look voom to construct dispersion
#a <- simulate_corr_vector(n_genes, 1, sigma)
simulate_corr_vector <- function(n_genes, 
                                 seed, sigma,
                                 log2FC = NULL, 
                                 diffexp_genes = NULL)
{

  # Generate random expression vector
  expression_vector <- MASS::mvrnorm(1, rep(10, n_genes), sigma)
  if (!is.null(diffexp_genes) & !is.null(log2FC)) {
    log2FC <- rep(0, n_genes)
    names(log2FC) <- gene_names
    
    # Taking genes from exact pathways and changing log2FC of them from 0 to 1
    log2FC[diffexp_genes] <- 1
    expression_vector <- expression_vector + log2FC
  }  
  return(expression_vector)
}

#for (n in c(4, 8, 10, 20)) {
#for (i in 1:5) {

pathways <- generated_pathways
pathways <- lapply(pathways, as.character)
n_genes <- length(gene_names)
n <- 5
i <- 1
print(i)
print(n)

diffexp_genes <- unique(c(pathways[["GOBP_ACTIN_FILAMENT_POLYMERIZATION"]]))


# Simulating clear vector of expressions 
matrix_expression_state_1 <- matrix(nrow = n_genes, ncol = n)
rownames(matrix_expression_state_1) <- gene_names

# simulate samples and append to the matrix
for (j in 1:n) {
  expression_state_1 <- simulate_corr_vector(n_genes,
                                             seed = i,
                                             sigma = sigma,
                                             )
  matrix_expression_state_1[, j] <- expression_state_1
}


matrix_expression_state_2 <- matrix(nrow = length(gene_names), ncol = n)
rownames(matrix_expression_state_2) <- gene_names
for (j in 1:n) {
  expression_state_2 <- simulate_corr_vector(n_genes,
                                             seed = i,
                                             sigma = sigma,
                                             log2FC = 1,
                                             diffexp_genes = diffexp_genes)
  matrix_expression_state_2[, j] <- expression_state_2
}

set.seed(i)


colnames(matrix_expression_state_1) <- rep('State_1', n) 
colnames(matrix_expression_state_2) <- rep('State_2', n) 
count_matrix <- cbind(matrix_expression_state_1, matrix_expression_state_2)

# Extracting pathways, which containing true differentially expressed genes(Which
# log2FC was changed)

pathways_containing_diffexp_genes <- names(pathways[lengths(lapply(pathways, intersect, diffexp_genes))])

# STARTING Limma
design <- data.frame(Group = factor(c(rep("group1", n), rep("group2", n))))
design_matrix <- model.matrix(~ Group, data = design)

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
result_table_fgsea_SE <- create_res(fgseaSERes, 'FGSEA_SE')
summary_table_SE <- rbind(summary_table_SE, result_table_fgsea_SE)

# Taking pathways, with pval < 0.05 after FGSEA_SE (function stored in report.R)
diffexpressed_pathways_fgsea_SE <- diffexpressed_pathways(result_table_fgsea_SE, 'padj')

# Calculating FP, TN, FN, TP and accuracy and adding it to summary report
report_fgsea_SE <- report('fgsea_SE', pathways_containing_diffexp_genes, 
                          diffexpressed_pathways_fgsea_SE, n, i, pathways)
summary_report <- rbind(report_fgsea_SE, summary_report)



# STARTING CAMERA

cameraRES <- camera_voom_m0(count_matrix, design_matrix, pathways)

result_table_camera <- create_res(cameraRES, 'camera')
summary_table_camera <- rbind(summary_table_camera, result_table_camera)


diffexpressed_pathways_camera <- diffexpressed_pathways(result_table_camera, 'padj')

report_camera <- report('CAMERA', pathways_containing_diffexp_genes,
                        diffexpressed_pathways_camera, n, i, pathways)
summary_report <- rbind(report_camera, summary_report)


# STARTING FORA


foraRes <- fora_res(de, pathways)


result_table_fora <- data.table(pathway = foraRes$pathway,
                                pval = foraRes$pval,
                                padj = p.adjust(foraRes$pval, method = 'BH'),
                                seed = i,
                                N = n,
                                Method = 'fora')
summary_table_fora <- rbind(summary_table_fora, result_table_fora)


diffexpressed_pathways_fora <- diffexpressed_pathways(result_table_fora, 'padj')

report_fora <- report('fora', pathways_containing_diffexp_genes, diffexpressed_pathways_fora,
                      n, i, pathways)
summary_report <- rbind(report_fora, summary_report)

# STARTING FGSEA MULTILEVEL  

fgseaRes <- fgseaMultilevel(pathways = pathways, 
                            stats)

# Creating result table for this N and seed and adding it to summary table 
result_table_fgsea <- create_res(fgseaRes, 'FGSEA')
summary_table_fgsea <- rbind(summary_table_fgsea, result_table_fgsea)

# Taking pathways, with pval < 0.05 after FGSEA_SE (function stored in report.R)
diffexpressed_pathways_fgsea <- diffexpressed_pathways(result_table_fgsea, 'padj')

# Calculating FP, TN, FN, TP and accuracy and adding it to summary report
report_fgsea <- report('fgsea', pathways_containing_diffexp_genes,
                       diffexpressed_pathways_fgsea, n, i, pathways)
summary_report <- rbind(report_fgsea, summary_report)

#}}
library(ggplot2)

summary_report$N <- factor(summary_report$N)


TPR <- ggplot(summary_report, aes(x=method, y=true_positive_rate, fill = N )) +  
  geom_boxplot() + theme_bw() + labs(title="True positive rate (1 pathway diffexpressed)")

TNR <- ggplot(summary_report, aes(x=method, y=true_negative_rate, fill = N )) +  
  geom_boxplot() + theme_bw() + labs(title="True negative rate (1 pathway diffexpressed)")

FPR <- ggplot(summary_report, aes(x=method, y=false_positive_rate, fill = N )) +  
  geom_boxplot() + theme_bw() + labs(title="False positive rate (1 pathway diffexpressed)")

FNR <- ggplot(summary_report, aes(x=method, y=false_negative_rate, fill = N )) +  
  geom_boxplot() + theme_bw() + labs(title="False negative rate (1 pathway diffexpressed)")

library("ggpubr")
figure <- ggarrange(TPR, TNR, FPR, FNR,
                    labels = c("A", "B", "C", 'D'),
                    ncol = 2, nrow = 2)
figure

#fraq <- ggplot(summary_report, aes(x = method, y = Fraction_of_differentially_expressed_pathways, fill = N)) +
#  geom_histogram() + theme_bw() + labs(title="True positive rate (1 pathway diffexpressed)")
#fraq

#fraq <- ggplot(summary_report, aes(x = method, fill = N)) +
#  geom_histogram(aes(y = Fraction_of_differentially_expressed_pathways)) +
#  theme_bw() +
#  labs(title="True positive rate (1 pathway diffexpressed)")
#fraq

FDR <- ggplot(summary_report, aes(x=method, y=false_discovery_rate, fill = N )) +  
  geom_boxplot() + theme_bw() + labs(title="False discovery (1 pathway diffexpressed)")

FDR

